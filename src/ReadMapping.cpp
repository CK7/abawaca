/*
 * ReadMapping.cpp
 *
 *  Created on: 31/Jan/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include "ReadMapping.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

namespace Bio {

/****************************************************************************************************************/
ReadMapping::ReadMapping(string _QNAME, size_t _FLAG, string _RNAME, size_t _POS, string CIGAR_indels, string _RNEXT, size_t _PNEXT, int _TLEN, Bio::DNAString _SEQ, string _QUAL, string CIGAR_mismatches) :
		QNAME(_QNAME), FLAG(_FLAG), RNAME(_RNAME), POS(_POS), RNEXT(_RNEXT), PNEXT(_PNEXT), TLEN(_TLEN), SEQ(_SEQ), QUAL(_QUAL)
{
	this->determine_snps(SEQ, QUAL, CIGAR_indels, CIGAR_mismatches, POS); 	
}

/****************************************************************************************************************/
ReadMapping::ReadMapping(const char* SAM_line)
{
        // HWI-ST330_0096:2:2:5723:2221#TTAGGC/2        16      NODE_107_length_1448229_cov_184.575043  816000  255     85M     *       0       0       
        // AAGAANAAGCAACAAAAGGTAATTTGACAGTTGCAGGTCATGATTTAATGAACATAAAAAATAAAGAAGTTCCTTATTTACGTCG
        // A<6BB#3917.-8?4CDD?:D86DCD6;EC=BCDA,D=D:C>A+CCAA;0B>+B7ECEE?BEEFECDCCD7>1:;EE8DCEF8CF        XA:i:0  MD:Z:5A46A32    NM:i:2

        // HWI-ST330_0096:2:2:8826:2326#TTAGGC/2        4       *       0       0       *       *       0       0
        // TTTTACAGGCACTAATGCACACAGAAAAGGTTTCCAGGCAATTTTCAATGCTTGCTATAAAGCATAAATACTTAATAAAGAATTATTTTCCTTCAGATTA
        // HHHHHHHHHHEHHHHHHHHHHHHHHHHHHHDGHFFEFFDGGHFHEFBFBEHHEBHGHHHHHBEHG?FGGEHBG7HGBFH?DC?FBAB?=AFEBBFA+9DG XM:i:0

	string CIGAR_indels, CIGAR_mismatches;
	vector<string>	fs;

	split('\t', SAM_line, fs);

	QNAME = fs[0];
	FLAG = atoi(fs[1].c_str());
	RNAME = fs[2];
	POS = atoi(fs[3].c_str())-1;
	CIGAR_indels = fs[5];
	if(!strcmp(fs[6].c_str(), "=")) {
		m_same_ref = true;
		RNEXT = RNAME;
	}
	else {
		m_same_ref = false;
		RNEXT = fs[6];
	}
	PNEXT = atoi(fs[7].c_str())-1;
	TLEN = atoi(fs[8].c_str());
	SEQ = Bio::DNAString(fs[9]);
	QUAL = fs[10];

	// Look for the mismatch CIGAR string 
	for(int i=11; i<fs.size(); i++) {
		size_t pos_mdz = fs[i].find("MD:Z:");
		if(pos_mdz != string::npos) {
			size_t e = fs[i].find(" ", pos_mdz);
			if(e==string::npos) {
				CIGAR_mismatches = fs[i].substr(pos_mdz);
			}
			else {
				CIGAR_mismatches = fs[i].substr(pos_mdz, (e-pos_mdz+1));
			}
			break;
		}
	}	

	this->determine_snps(SEQ, QUAL, CIGAR_indels, CIGAR_mismatches, POS); 
}

/****************************************************************************************************************/
// Addition on 29/12/13 (to be compatible with bowtie2):
// CIGAR_indels is field #5 in the sam format
// CIGAR_mismatches is the field that starts with MD:Z.
void ReadMapping::determine_snps(const DNAString& dna_str, const string& quality_str, string CIGAR_indels, string CIGAR_mismatches, size_t ref_pos)
{
	if(CIGAR_mismatches.size() == 0)
		return;

	// This regexp was taken from the SAM manual
	// /^MD:Z:([0-9]+)(([ACGTN]|\^[ACGTN]+)[0-9]+)*$/)
	string::const_iterator p = CIGAR_mismatches.begin();
	if((*p!='M') || (*(p+1)!='D') || (*(p+2)!=':') || (*(p+3)!='Z') || (*(p+4)!=':')) {
		throw Illegal_mapping(CIGAR_mismatches, "Illegal SNPs description, expected MD:Z: at the beginning of CIGAR_mismatches format but found something different", __PRETTY_FUNCTION__);
	}
	p += 5;

	// Read the offset to the first SNP
	size_t	read_pos = 0, ref_deleted = 0;
	while((*p>='0') && (*p<='9')) {
		read_pos *= 10;
		read_pos += (*(p++)-'0');
	}
	// Set the offset to 0 
	read_pos--;
	// Keep on reading the next snps until we reach the end
	while(p != CIGAR_mismatches.end()) {
		// First we expect the ref's char
		if((*p!='^') && ((*p<'A') || (*p>'Z')))
			throw Illegal_mapping(CIGAR_mismatches, string("Illegal SNPs description (2), expected CIGAR_mismatches format but found something different"), __PRETTY_FUNCTION__);

		if(*p=='^') {
			p++;
			while(((*p<'0') || (*p>'9')) && (p != CIGAR_mismatches.end())) {
				ref_deleted++;
				m_snps.push_back(SNP(SNP::DELETION, ref_pos+read_pos+ref_deleted, read_pos, *p, '-', 40));
				p++;
			}
		}
		else {
			read_pos++;
			m_snps.push_back(SNP(SNP::MISMATCH, ref_pos+read_pos+ref_deleted, read_pos, *p, dna_str.get(read_pos), quality_str.at(read_pos)));
			p++;
		}

		// There must be at least one more number
		if((*p < '0') || (*p > '9'))
			throw Illegal_mapping(CIGAR_mismatches, "Illegal SNPs description, expeected CIGAR_mismatches format but found something different", __PRETTY_FUNCTION__);
		unsigned int o = 0;
		while((*p>='0') && (*p<='9')) {
			o *= 10;
			o += (*(p++)-'0');
		}
		read_pos += o;
	}

	// Now we have to go through the insertion/deletion CIGAR string and correct mismatches accordingly
	// Examples: 6M1D9M1I132M, 31M1I106M	
	p = CIGAR_indels.begin();
	// Make sure that the end of the string is legal (so we don't have to check it all the time in the loop below)
	char ch = *(CIGAR_indels.rbegin());
	if((ch!='I') && (ch!='D') && (ch!='M')) {
		throw Illegal_mapping(CIGAR_indels, string("Illegal CIGAR-indels description (3)"), __PRETTY_FUNCTION__);
	}

	int read_pos2 = 0;
	size_t ref_pos2 = ref_pos;
	while(p!=CIGAR_indels.end()) {
		size_t o = 0;
		while((*p >= '0') && (*p <= '9')) {
			o *= 10;
			o += (*(p++)-'0');
		}
		if(o == 0) {
			throw Illegal_mapping(CIGAR_indels, string("Illegal CIGAR-indels description (4)"), __PRETTY_FUNCTION__);
		}
		vector<SNP>::iterator ps;
		size_t pos;
		switch (*(p++)) {
			case 'M':
				read_pos2 += o;
				ref_pos2 += o;
				break;
			case 'D':
				// Nothig to do here - this case was already taken care of in CIGAR_mismatches
				ref_pos2 += o;
				break;
			case 'I':
				for(ps=m_snps.begin(); ps!=m_snps.end(); ps++) {
					if(ps->read_pos >= read_pos2) {
						ps->read_pos += o;
					}
				}
				read_pos += o;

				for(ps=m_snps.begin(); (ps!=m_snps.end()) && (ps->read_pos < read_pos2); ps++)
					;
				// At this point we don't provide any information about the insetion except that it exists
				m_snps.insert(ps, SNP(SNP::INSERTION, ref_pos2-1, read_pos2-1, '-', '-', '-'));
				read_pos2 += o;

				break;
			default:
				throw Illegal_mapping(CIGAR_indels, string("Unexpected CIGAR-indels operation found"), __PRETTY_FUNCTION__);
		}
	}

	// Sanity check
	if(((read_pos+1) != dna_str.size()) || (read_pos2 != dna_str.size()))
		throw Illegal_mapping(CIGAR_indels+string(", ")+CIGAR_mismatches, "snps and offsets do not sum to expected read length", __PRETTY_FUNCTION__);

}

}
