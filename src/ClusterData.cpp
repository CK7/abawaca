#include "ClusterData.h"
#include <iostream>
#include <exception>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClusterData::ClusterData(const ClusterData& cluster_data, const unordered_set<size_t>& dp_subset) : dps(dp_subset)
{
	dimensions.push_back(NULL);
        for(size_t i=1; i<=cluster_data.ndimensions(); i++)  {
                dimensions.push_back(new Dimension(cluster_data.get_dimension(i), dps));
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClusterData::~ClusterData()
{
	for(vector<Dimension*>::iterator vit=dimensions.begin(); vit!=dimensions.end(); vit++)
		if(*vit != NULL)
			delete(*vit);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClusterData::ClusterData(string esom_lrn_file, const ScafDpData& scaf_db)
{
        FILE* fp = fopen(esom_lrn_file.c_str(), "r");

        if(!fp) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: could not read " << esom_lrn_file << endl << endl;
                exit(-1);
        }

        char line[8192] = {};

	/********** Read # of dps **********/
	// % 8430
        if(!fgets(line, 8192, fp)) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: " << esom_lrn_file << " is empty" << endl << endl;
                exit(-1);
        }

        // % 8430
        unsigned int ndps;
        char    c;
        if((sscanf(line, "%c %u", &c, &ndps) != 2) || (c != '%')) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: unexpected line 1 in " << esom_lrn_file << ": " << endl << line << endl << endl;
                exit(-1);
        }

	/********** Read # of dimensions **********/
	// % 12
        if(!fgets(line, 8192, fp)) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: failed to read line 2 in " << esom_lrn_file << endl << endl;
                exit(-1);
        }

        unsigned int ndimensions;
        if((sscanf(line, "%c %u", &c, &ndimensions) != 2) || (c != '%')) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: unexpected line 2 in " << esom_lrn_file << ": " << endl << line << endl << endl;
                exit(-1);
        }
	ndimensions--;		// Count in the file contains the dp column

	/********** Read weights for dimensions. We don't use it but just verify that the number found is right **********/
	// % 9	1	1	1	1	1	1	1	1	1	1	1
        if(!fgets(line, 8192, fp)) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: failed to read line 3 in " << esom_lrn_file << endl << endl;
                exit(-1);
        }

        unsigned int t = 1;
	while(strrchr(line, '\t') != NULL) {
		t++;
		*strrchr(line, '\t') = 0;
	}
        if(t != (ndimensions+1)) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: unexpected line 3 in " << esom_lrn_file << ": number of fileds is differemt that " << ndimensions << endl << line << endl << endl;
                exit(-1);
        }

	/********** Read dimension names **********/
	// % key	74-1	74-10	74-11	74-12	74-2	74-3	74-4	74-5	74-6	74-7	74-9
        if(!fgets(line, 8192, fp)) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: failed to read line 4 in " << esom_lrn_file << endl << endl;
                exit(-1);
        }

	// Chomp
	while((strlen(line)>0) && isspace(*(line+strlen(line)-1)))
		*(line+strlen(line)-1) = 0;

	dimensions.push_back(NULL);
	// We skip the first column which contains the dp title
	char *p1=strchr(line, '\t')+1, *p2=strchr(p1, '\t');
	size_t	i = 1;
	while(p2) {
		*p2 = 0;
		dimensions.push_back(new Dimension(i++, string(p1)));
		p1 = p2+1;
		p2=strchr(p1, '\t');
	}
	// Last dimension is not accompanied by a tab
	dimensions.push_back(new Dimension(i, string(p1)));

	if(i != ndimensions) {
                cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: number of dimensions found in line 4 (" << i << ") is different than expected (" << ndimensions << ") in file " << esom_lrn_file << ":" << endl << line << endl << endl;
                exit(-1);
	}

	/********** And finally, read the values themselves **********/	
	size_t entry=0;
	map<size_t, string>	data;

	while(fgets(line, 8192, fp)) {
		// Chomp
		while((strlen(line)>0) && isspace(*(line+strlen(line)-1)))
			*(line+strlen(line)-1) = 0;

		// Skip empty or comment lines
		if((line[0] == '%') || (strlen(line) == 0))
			continue;

		// The following means that we have more lines than expected by the % dps information
		if(entry >= ndps) {
			cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: number of datapoint lines found is higher than expected (" << ndps << "), file " << esom_lrn_file << endl << endl;
			exit(-1);
		}

		char *p1, *p2;
		size_t t=1;
		// Count the number of values and make sure it's what we would expect
		for(p1=line, p2=strchr(p1, '\t'); p2!=NULL; p1=p2+1, p2=strchr(p1, '\t'))
			t++;
		if(t != (ndimensions+1)) {
			cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: number of dimensions found is " << t << ", expected to find " << ndimensions << ": " << line << endl << endl; 
			exit(-1);
		}

		// 1	0.181737916437441	0	0.00990659188488923	0.551726686306843	0.0250488966425045	0.0706148967876072	0	0	0.00324475022645163	0.0147897529305731	0.142930508783691
		p1=line;
		p2=strchr(p1, '\t');
		*p2=0;
		int dp = scaf_db.dp_name2id(atoi(p1));

		if(dp != 0) {
			data.insert(pair<size_t, string>(dp, string(p2+1)));
		}
		entry++;
	}

	for(auto it=data.begin(); it!=data.end(); it++) {
		size_t dp = it->first;
		dps.insert(dp);
		int p1, p2, t;
		for(t=1, p1=0, p2=it->second.find('\t', 0); p2!=string::npos; t++, p1=p2+1, p2=p2=it->second.find('\t', p1)) {
			dimensions[t]->insert(dp, atof(it->second.substr(p1, p2-p1).c_str()));
		}
		// Last value
		dimensions[t]->insert(dp, atof(it->second.substr(p1).c_str()));
	}
	if(entry != ndps) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: number of datapoint lines found (" << entry << ") is different than expected (" << ndps << "), file " << esom_lrn_file << endl << endl;
		exit(-1);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const Dimension& ClusterData::get_dimension(string dname) const
{
	auto it=dimensions.begin();
	for(++it; it!=dimensions.end() && ((*it)->get_dimension_name() != dname); it++)
		;

	if(it == dimensions.end()) {
		char s[2048];
		sprintf(s, "%s (%d): could not find dimension %s", __FILE__, __LINE__, dname.c_str());
		throw invalid_argument(s);
	}

	return **it;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const Dimension& ClusterData::get_dimension(size_t index) const
{
	if((index == 0) || (index >= dimensions.size())) {
		char s[2048];
		sprintf(s, "%s (%d): illegal dimension index %lu, should be in the range of [0, %lu]", __FILE__, __LINE__, index, dimensions.size()-1);
		throw invalid_argument(s);
	}

	return *(dimensions[index]);
}
