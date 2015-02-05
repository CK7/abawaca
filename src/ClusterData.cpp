#include "ClusterData.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClusterData::ClusterData(const ClusterData& cluster_data, const set<size_t>& dp_subset) : dps(dp_subset)
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
ClusterData::ClusterData(const char* esom_lrn_file)
{
        FILE* fp = fopen(esom_lrn_file, "r");

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
		int dp = atoi(p1);
		dps.insert(dp);
		for(t=1, p1=p2+1, p2=strchr(p1, '\t'); p2!=NULL; t++, p1=p2+1, p2=strchr(p1, '\t')) {
			*p2 = 0;
			dimensions[t]->insert(dp, atof(p1));
		}
		// Last value
		dimensions[t]->insert(dp, atof(p1));
		entry++;
	}
	if(entry != ndps) {
		cerr << __FILE__ << " (" << __LINE__ << "): Fatal error: number of datapoint lines found is different than expected (" << ndps << "), file " << esom_lrn_file << endl << endl;
		exit(-1);
	}
}
