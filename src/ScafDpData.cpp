#include <iostream>
#include <stdio.h>
#include <string.h>
#include <exception>
#include "ScafDpData.h"
#include "SeqIORead_fasta.h"

using namespace Bio;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ScafDpData::ScafDpData(string esom_names_file, string fasta_file, string info_file)
{
	FILE* fp = NULL;
	char line[8192] = {};

	if(!(fp = fopen(info_file.c_str(), "r"))) {
		throw invalid_argument(string(__FILE__) + " (" + to_string(__LINE__) + "): could not read " + info_file);
	}

	map<string, double>	scaf_name2cvg;
	while(fgets(line, 8192, fp)) {
		// gwd1_scaffold_0 556922  17.599  0.363   0
		char seq_name[512];
		double cvg, gc;
		size_t Ns, length;
		if(sscanf(line, "%s %lu %lf %lf %lu", seq_name, &length, &cvg, &gc, &Ns) != 5) {
			char s[2048];
			sprintf(s, "%s (%d): Illegal line in file %s:\n%s", __FILE__, __LINE__, info_file.c_str(), line);
			throw invalid_argument(s);
		}
		scaf_name2cvg.insert(pair<string, double>(seq_name, cvg));
	}
	fclose(fp);

	if(!(fp = fopen(esom_names_file.c_str(), "r"))) {
		char s[2048];
		sprintf(s, "%s (%d): could not read %s", __FILE__, __LINE__, esom_names_file.c_str());
		throw invalid_argument(s);
	}

	if(!fgets(line, 8192, fp)) {
		char s[2048];
		sprintf(s, "%s (%d): file %s is empty", __FILE__, __LINE__, esom_names_file.c_str());
		throw invalid_argument(s);
	}

	// % 8430
	unsigned int t1, t2;
	char	c;
	if((sscanf(line, "%c %u", &c, &t1) != 2) || (c != '%')) {
		char s[2048];
		sprintf(s, "%s (%d): unexpected line in file %s:\n%s", __FILE__, __LINE__, esom_names_file.c_str(), line);
		throw invalid_argument(s);
	}

	map<string, set<DataPoint> >	scaf2dps;
	while(fgets(line, 8192, fp)) {
		// 1    gwd1_scaffold_20283_1   gwd1_scaffold_20283:(1,2072), 1 segment(s), 2072/2072 non-N bps
		// 1	gwd1_scaffold_0_1	gwd1_scaffold_0(1, 2003), 2003/2003 non-Ns bps
		size_t 	dp=0, start=0, end=0, non_Ns=0, length=0, nsegments=0;
		char	scaf[512]{}, gene_name[512]={}, *p1, *p2, *p3, *p4, *p5;

		p1 = strchr(line, ':');
		p2 = strchr(line, '(');
		p3 = strchr(line, ')');
		p4 = strchr(line, '/');
		p5 = strchr(line, ',');
		if(!p1 || !p2 || !p3 || !p4 || !p5) {
			char s[2048];
			sprintf(s, "\n%s (%d): unexpected line in file %s:\n%s\n", __FILE__, __LINE__, esom_names_file.c_str(), line);
			throw invalid_argument(s);
		}
		*p1 = *p2 = *p3 = *p4 = *p5 = ' ';
		if((sscanf(line, "%lu %s %s %lu %lu , %lu segment(s), %lu %lu non-N bps", &dp, gene_name, scaf, &start, &end, &nsegments, &non_Ns, &length) != 8) &&
		   (sscanf(line, "%lu %s %s %lu %lu , %lu %lu non-N bps", &dp, gene_name, scaf, &start, &end, &non_Ns, &length) != 7)) {
			char s[2048];
			sprintf(s, "\n%s (%d): unexpected line in file %s:\n%sdp=%lu, gene_name=%s, scaf=%s, start=%lu, end=%lu, nsegments=%lu, non_Ns=%lu, length=%lu\nMake sure that your esom files were generated using prepare_esom_files.pl v1.05 or later\n",
				 __FILE__, __LINE__, esom_names_file.c_str(), line, dp, gene_name, scaf, start, end, nsegments, non_Ns, length);
			throw invalid_argument(s);
		}
		auto sit = scaf2dps.find(scaf);
		if(sit == scaf2dps.end()) {
			scaf2dps.insert(pair<string, set<DataPoint> >(scaf, set<DataPoint>()));
			sit = scaf2dps.find(scaf);
		}
		sit->second.insert(DataPoint(0, dp, 0, start, end, length-non_Ns));
	}
	fclose(fp);

	dps.resize(1);
	for(auto sit = scaf2dps.begin(); sit != scaf2dps.end(); sit++) {
		if(sit->second.size() == 1)
			continue;
		scaf_name2id_map.insert(pair<string, size_t>(sit->first, scaf_name2id_map.size()+1));
		for(auto dpit=sit->second.begin(); dpit != sit->second.end(); dpit++) {
			dps.push_back(DataPoint(dps.size(), dpit->get_name(), scaf_name2id_map.size(), dpit->get_start(), dpit->get_end(), dpit->get_Ns()));
			dp_name2id_map.insert(pair<size_t, size_t>(dps.back().get_name(), dps.back().get_iid()));
		}
	}

	scafs.resize(scaf_name2id_map.size()+1);
	Bio::SeqIORead_fasta<Bio::DNASequence>	reader(fasta_file);
	Bio::DNASequence* seq_obj;

	while((seq_obj = reader.next_seq()) != NULL) {
		auto it=scaf_name2id_map.find(seq_obj->display_id());
		auto mit = scaf_name2cvg.find(seq_obj->display_id());
		double cvg = (mit == scaf_name2cvg.end())? 0 : mit->second;
		if(it != scaf_name2id_map.end()) {
			scafs[it->second] = Seq(it->second, seq_obj, cvg);
		}
		else {
			delete(seq_obj);
		}
	}

	auto it=dps.begin();
	for(++it; it!=dps.end(); it++)
		scafs[it->get_scaf_id()] += *it;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t ScafDpData::dp2scaf(size_t dp) const
{
	if((dp == 0) || (dp >= dps.size())) {
		char s[2048];
		sprintf(s, "%s (%d): Invalid dp %lu", __FILE__, __LINE__, dp);
		throw invalid_argument(s);
	}
	return dps[dp].get_scaf_id();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t ScafDpData::ndps(size_t scaf_id) const
{
	if((scaf_id == 0) || (scaf_id >= scafs.size())) {
		char s[2048];
		sprintf(s, "%s (%d): Invalid scaf_id %lu", __FILE__, __LINE__, scaf_id);
		throw invalid_argument(s);
	}
        return scafs[scaf_id].get_scaf_dps().size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string ScafDpData::scaf_id2name(size_t scaf_id) const
{
	if((scaf_id == 0) || (scaf_id >= scafs.size())) {
		char s[2048];
		sprintf(s, "%s (%d): Invalid scaf_id %lu", __FILE__, __LINE__, scaf_id);
		throw invalid_argument(s);
	}
        return scafs[scaf_id].get_seq_obj()->display_id();
}
