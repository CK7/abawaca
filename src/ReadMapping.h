/*
 * ReadMapping.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#ifndef READMAPPING_H_
#define READMAPPING_H_

#include "common.h"
#include "String.h"
#include <memory>
#include <iostream>
namespace Bio {

/***************************************************************************************************************************
 * ReadMapping
 * Stores information about the mapping of a read to a reference sequence.
 * Uses 0 offset (like all other classes).
 * Class follows specification for the SAM format, see here: http://samtools.github.io/hts-specs/SAMv1.pdf
 ***************************************************************************************************************************/
class ReadMapping {
public:
	class Illegal_mapping;
public:
	struct SNP;
	typedef enum {forward, reverse} Direction;
public:
	// ref_pos is expected to be a 0-offset position.
	ReadMapping(const char* SAM_line);
 	ReadMapping(string _QNAME, size_t _FLAG, string _RNAME, size_t _POS, string CIGAR_indels, string _RNEXT, size_t _PNEXT, int _TLEN, Bio::DNAString _SEQ, string _QUAL, string CIGAR_mismatches);
	virtual 		~ReadMapping()					{/*std::cout << "Destroying " << RNAME << std::endl;*/}
	const string&		ref_name() const				{return RNAME;}
	const string&		pair_ref_name() const				{return RNEXT;}
	bool			pair_mapped_to_same_ref() const			{return (this->mapped() && this->pair_mapped() && m_same_ref);}
	const string&		read_name() const				{return QNAME;}
	const size_t		read_length() const				{return SEQ.size();}
	const Bio::DNAString	mapped_read() const				{return SEQ;}
	const string&		quality_string() const				{return QUAL;}
	size_t			ref_pos() const					{return POS;}
	size_t			pair_ref_pos() const				{return PNEXT;}
	int			insert_size() const				{return TLEN;}
	Direction		direction() const				{return (FLAG & 0x10)? reverse : forward;}
	Direction		pair_direction() const				{return (FLAG & 0x20)? reverse : forward;}
	size_t			num_snps() const				{return m_snps.size();}
	// Index may be 0 .. (num_snps()-1)
	const SNP&		snp(size_t index) const				{return m_snps[index];}
	bool			unmapped() const				{return (FLAG & 0x4);}
	bool			pair_unmapped() const				{return (FLAG & 0x8);}
	bool			mapped() const					{return !(FLAG & 0x4);}
	bool			pair_mapped() const				{return !(FLAG & 0x8);}
	bool			multiple_hits() const				{return (FLAG & 0x100); }
	bool			first_in_segment() const			{return (FLAG & 0x40);}
	bool			last_in_segment() const				{return (FLAG & 0x80);}
protected:
	void 			determine_snps(const DNAString& dna_str, const string& quality_str, string CIGAR_indels, string CIGAR_mismatches, size_t ref_pos);
protected:
	string			QNAME;
	size_t			FLAG;
	string			RNAME;
	size_t			POS;
	string			RNEXT;
	size_t			PNEXT;
	int			TLEN;
	Bio::DNAString		SEQ;
	string			QUAL;
	bool			m_same_ref;
	vector<SNP>		m_snps;
};

//typedef std::tr1::shared_ptr<ReadMapping> ReadMappingPtr;
typedef std::shared_ptr<ReadMapping> ReadMappingPtr;

/****************************************************************************************************************/
struct ReadMapping::SNP {
	typedef enum {INSERTION=0x01, DELETION=0x02, MISMATCH=0x04} SNP_type;
	SNP_type		snp_type;
	size_t			ref_pos, read_pos;
	char			ref_char, read_char;
	unsigned long		read_quality;
	SNP(SNP_type _snp_type, size_t _ref_pos, size_t _read_pos, char _ref_char, char _read_char, unsigned long _read_quality) :
		snp_type(_snp_type), ref_pos(_ref_pos), read_pos(_read_pos), ref_char(_ref_char), read_char(_read_char), read_quality(_read_quality)
		{}
	SNP() : ref_pos(0), read_pos(0), ref_char('-'), read_char('-'), read_quality(0)
		{}
};

/****************************************************************************************************************/
class ReadMapping::Illegal_mapping : public exception {
	string msg;
public:
	Illegal_mapping(string mapping_str, string err_msg, const char* func)
		{
			msg = string("Fatal error:\t") + err_msg + string("\nMapping string:\t") + mapping_str + string("\n");
			msg += string("Function:  \t") + string(func) + string("\n");

		}
	~Illegal_mapping() throw()
		{}
	virtual const char* what() const throw()
		{return msg.c_str();}
};

}
#endif /* READMAPPING_H_ */
