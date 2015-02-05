#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ScafDpData.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ScafDpData::ScafDpData(const char* esom_names_file) : dp2scaf_array(NULL), scaf_id2ndps_array(NULL), num_scafs(0), num_dps(0)
{
	FILE* fp = fopen(esom_names_file, "r");

	if(!fp) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: could not read " << esom_names_file << endl << endl;
		exit(-1);
	}

	char line[8192] = {};
	if(!fgets(line, 8192, fp)) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: " << esom_names_file << " is empty" << endl << endl;
		exit(-1);
	}

	// % 8430
	unsigned int t1, t2;
	char	c;
	if((sscanf(line, "%c %u", &c, &t1) != 2) || (c != '%')) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: unexpected line in " << esom_names_file << ": " << endl << line << endl << endl;
		exit(-1);
	}
	num_dps = t1;

	if(!(dp2scaf_array = (size_t*)calloc(num_dps+1, sizeof(size_t)))) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: could not allocate space for " << num_dps << " datapoints as specified in " << esom_names_file << endl << endl;
		exit(-1);
	}

	while(fgets(line, 8192, fp)) {
		// 1	311210.1119.1809_1	311210.1119.1809
		size_t 	dp;
		char	scaf[512], desc[512];
		if((sscanf(line, "%u %s %s", &t1, scaf, desc) != 3) || (strrchr(scaf, '_') == NULL)) {
			cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: unexpected line in " << esom_names_file << ": " << endl << line << endl << endl;
			exit(-1);
		}
		dp = t1;
		char* p;
		for(p=scaf+strlen(scaf)-1; *p!='_'; p--) {
			if((*p<'0') || (*p>'9')) {
				cerr << __FILE__ << " (" << __LINE__ << "): Unexpected postfix in datapoint name (second field) of this line, should be digits only: " << line << endl;
				exit(-1);
			}
		}
		*p = 0;

		map<string, size_t>::const_iterator it = scaf_name2id_map.find(scaf);
		if(it == scaf_name2id_map.end()) {
			scaf_name2id_map[scaf] = ++num_scafs;
			it = scaf_name2id_map.find(scaf);
		}
		dp2scaf_array[dp] = it->second;
	}
	fclose(fp);

	if(!(scaf_id2ndps_array = (size_t*)calloc(num_scafs+1, sizeof(size_t)))) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: could not allocate space for " << num_dps << " datapoints as specified in " << esom_names_file << endl << endl;
		exit(-1);
	}

	scaf_id2name_array.resize(num_scafs+1);
	for(map<string, size_t>::const_iterator it = scaf_name2id_map.begin(); it!=scaf_name2id_map.end(); it++)
		scaf_id2name_array[it->second] = it->first;


	for(size_t dp=1; dp<=num_dps; dp++)
		scaf_id2ndps_array[dp2scaf_array[dp]]++;

}
