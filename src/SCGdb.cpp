#include "SCGdb.h"
#include <exception>
#include <string.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t	SCGdb::num_unique_scgs(const set<size_t>& scafs) const
{
	set<string> scgs_found;

	for(auto it=scafs.begin(); it!=scafs.end(); it++) {
		auto mit = scaf2scgs.find(*it);
		if(mit == scaf2scgs.end())
			continue;
		scgs_found.insert(mit->second.begin(), mit->second.end());
	}

	return scgs_found.size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SCGdb::get_scg2count(const set<size_t>& scafs, map<string, size_t>& result) const
{
	result.clear();

	for(auto it=scafs.begin(); it!=scafs.end(); it++) {
		auto mit = scaf2scgs.find(*it);
		if(mit == scaf2scgs.end())
			continue;
		for(auto scgit = mit->second.begin(); scgit!=mit->second.end(); scgit++) {
			auto rit = result.find(*scgit);
			if(rit == result.end()) {
				result.insert(pair<string, size_t>(*scgit, 0));
				rit = result.find(*scgit);
			}
			rit->second++;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double	SCGdb::average_num_copies_for_unique_scgs(const set<size_t>& scafs) const
{
	set<string> scgs_found;
	size_t n = 0;

	for(auto it=scafs.begin(); it!=scafs.end(); it++) {
		auto mit = scaf2scgs.find(*it);
		if(mit == scaf2scgs.end())
			continue;
		scgs_found.insert(mit->second.begin(), mit->second.end());
		n += mit->second.size();
	}

	return (n == 0)? 0 : double(n)/scgs_found.size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SCGdb::load_scgs()
{
	FILE* fin = fopen(scg_list_file.c_str(), "r");
	if(fin == NULL) {
		char	s[2048];
		sprintf(s, "Failed to read file %s", scg_list_file.c_str());
		throw invalid_argument(s);
	}

	char	scg[1024], gene[1024];

	total_num_scgs = 0;
	while(fscanf(fin, "%s", scg) == 1) {
		total_num_scgs++;
	}
	fclose(fin);

	// An empty SCG database
	if(gene2scg_file.size() == 0)
		return;

	fin = fopen(gene2scg_file.c_str(), "r");
	if(fin == NULL) {
		char	s[2048];
		sprintf(s, "Failed to read file %s", gene2scg_file.c_str());
		throw invalid_argument(s);
	}

	while(fscanf(fin, "%s %s", gene, scg) == 2) {
		char* p = strrchr(gene, '_');
		if(p == NULL) {
			char	s[2048];
			sprintf(s, "Unexpected gene name %s, expected prodigal format: <name>_<integer>\nfile: %s", gene, scg_list_file.c_str());
			throw invalid_argument(s);
		}

		*p = 0;
		size_t scaf = scaf_db.scaf_name2id(gene);
		// The following can happen if a sequence is in the file that was used for the SCG analysis but was not in the file that was binned.
		// We are fine with that.
		if(scaf == 0) {
			continue;
		}

		int gene_id = atoi(p+1);

		if(gene_id == 0) {
			string scaf_name = gene;
			*p = '_';
			char	s[2048];
			sprintf(s, "Could not identify gene id in %s (recovered scaf=%s (%lu), gene_id=%s). Expecting prodigal format gene names: <name>_<integer>\nfile: %s",
				gene, scaf_name.c_str(), scaf, p+1, scg_list_file.c_str());
			throw invalid_argument(s);
		}
		auto it = scaf2scgs.find(scaf);
		if(it == scaf2scgs.end()) {
			scaf2scgs.insert(pair<size_t, set<string> >(scaf, set<string>()));
			it = scaf2scgs.find(scaf);
		}
		it->second.insert(scg);
	}
	fclose(fin);
}
