#ifndef SCAFDPDATA_H
#define	SCAFDPDATA_H

#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ScafDpData {
public:
	ScafDpData(const char* esom_names_file);
	~ScafDpData()					{free(dp2scaf_array); free(scaf_id2ndps_array);}
	size_t	nscafs() const				{return num_scafs;}
	size_t	ndps() const				{return num_dps;}
	size_t	ndps(size_t scaf_id) const		{return scaf_id2ndps_array[scaf_id];}
	size_t	dp2scaf(size_t dp) const		{return dp2scaf_array[dp];}
	string	scaf_id2name(size_t scaf) const		{return scaf_id2name_array[scaf];}
	size_t	scaf_name2id(string scaf) const		{if(scaf_name2id_map.find(scaf)!=scaf_name2id_map.end()) return scaf_name2id_map.find(scaf)->second; return 0;}
protected:
	size_t*			dp2scaf_array;
	size_t*			scaf_id2ndps_array;
	vector<string>		scaf_id2name_array;
	map<string, size_t>	scaf_name2id_map;
	size_t			num_scafs;
	size_t			num_dps;
};

#endif	// SCAFDPDATA_H
