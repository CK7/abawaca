#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "misc.h"
#include "ScafDpData.h"
#include "ClusterData.h"
#include "ClusterSeparator.h"
#include "ClusterWriter.h"

using namespace std;

size_t get_fasta(const set<size_t>& scafs, const ScafDpData& scaf_db, const char* fasta_out, const char* fasta_file);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cluster* init_cluster(const ClusterData& cluster_db, const ScafDpData& scaf_db)
{
	Cluster* cluster = new Cluster(scaf_db);
	for(size_t i=1; i<=scaf_db.ndps(); i++)
		*cluster += i;
	for(size_t i=1; i<=scaf_db.nscafs(); i++)
		cluster->add_assigned_scaf(i);

	return cluster;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void usage(const char* prog_name)
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char** argv) 
{
	if(argc < 5) {
		fprintf(stderr, "\nUsage: %s <.names file> <.lrn file> <fasta-file> <out-dir>\n\n", argv[0]);
		exit(-1);
	}

	if(argc < 5) {
		fprintf(stderr, "\nUsage: %s -n <.names file> -l <.lrn file> -f <fasta-file> -o <out-dir>\n\n", argv[0]);
		exit(-1);
	}



	ScafDpData			scaf_db(argv[1]);
	ClusterData			clustering_db(argv[2]);
	map<size_t, Cluster*>		clusters;
	size_t				total_iterations = 1;
	map<size_t, size_t>		dp2cluster;
	map<size_t, size_t>		scaf2cluster;
	const char* 			root_dir = argv[4];
	const char*			fasta_file = argv[3];
	char				cluster_dir[1024];
	char				fasta_dir[1024];
	ofstream			log;
	ofstream			summary;
	ofstream			ofs_dp2cluster;
	ofstream			ofs_scaf2cluster;
	char 				line[1024] = {};

	sprintf(cluster_dir, "%s/clusters", root_dir);
	sprintf(fasta_dir, "%s/final-clusters", root_dir);
	sprintf(line, "mkdir %s; mkdir %s; mkdir %s", root_dir, cluster_dir, fasta_dir);
	system(line);
	sprintf(line, "%s/log", root_dir);
	log.open(line);
	sprintf(line, "%s/summary.txt", root_dir);
	summary.open(line);
	sprintf(line, "%s/dp2cluster.txt", root_dir);
	ofs_dp2cluster.open(line);
	sprintf(line, "%s/scaf2cluster.txt", root_dir);
	ofs_scaf2cluster.open(line);


	summary << "Cluster\t# scafs\t# dps\tTotal bps" << endl;

	clusters[total_iterations++] = init_cluster(ClusterData(argv[2]), scaf_db);

	// We will keep iterating over all items in clustering_dbs until it is empty.
	while(clusters.size() > 0) {
		map<size_t, Cluster*>::iterator 	mit = clusters.begin();
		Cluster* 				curr_cluster = mit->second; 
		ClusterData&    			curr_clustering_db = *(new ClusterData(clustering_db, curr_cluster->get_dps()));
		size_t 					iteration = mit->first;

		clusters.erase(mit);

		cerr << endl << '[' << get_time() << ']' << " Checking cluster " << iteration << " (" << curr_cluster->ndps() << " datapoints)" << endl;
		log << endl << '[' << get_time() << ']' << " Checking cluster " << iteration << " (" << curr_cluster->ndps() << " datapoints)" << endl;

		ClusterSeparatorBySensitivitySpecificity separator(scaf_db, curr_clustering_db);

		separator.separate();

		if(separator.get_cluster1() == NULL) {
			cerr << '[' << get_time() << ']' << " Cluster " << iteration << " is a terminal cluster" << endl; 
			log << '[' << get_time() << ']' << " Cluster " << iteration << " is a terminal cluster" << endl; 
			for(set<size_t>::const_iterator sit = curr_cluster->dps_begin(); sit != curr_cluster->dps_end(); sit++)
				dp2cluster[*sit] = iteration;
			for(set<size_t>::const_iterator sit = curr_cluster->assigned_scafs_begin(); sit != curr_cluster->assigned_scafs_end(); sit++)
				scaf2cluster[*sit] = iteration;
			char fasta_out[1024];

			sprintf(fasta_out, "%s/%lu.fasta", fasta_dir, iteration);
			size_t total_bps = get_fasta(curr_cluster->get_assigned_scafs(), scaf_db, fasta_out, fasta_file);
			summary << iteration << '\t' << curr_cluster->get_assigned_scafs().size() << '\t' << curr_cluster->ndps() << '\t' << total_bps << endl; 
		}
		else {
			cerr << '[' << get_time() << ']' << " Best separation:\t" 
				<< "dimension " << separator.get_separating_dimension() << " (" << curr_clustering_db.get_dimension(separator.get_separating_dimension()).get_dimension_name() << ")\t"
				<< "Value " << separator.get_separating_value() << "\t"
				<< "Specificity " << separator.get_best_specificity() << "\t"
				<< "Sensitivity " << separator.get_best_sensitivity() << endl;
			cerr << '[' << get_time() << ']' << " Cluster " << iteration << " will be separated to the following two clusters" << endl; 
			log << '[' << get_time() << ']' << " Best separation:\t" 
				<< "dimension " << separator.get_separating_dimension() << " (" << curr_clustering_db.get_dimension(separator.get_separating_dimension()).get_dimension_name() << ")\t"
				<< "Value " << separator.get_separating_value() << "\t"
				<< "Specificity " << separator.get_best_specificity() << "\t"
				<< "Sensitivity " << separator.get_best_sensitivity() << endl;
			log << '[' << get_time() << ']' << " Cluster " << iteration << " will be separated to the following two clusters" << endl; 

			ClusterEsomWriter writer(cluster_dir, scaf_db);
			char new_cluster_name[128];

			sprintf(new_cluster_name, "%lu", total_iterations);
			ClusterData cluster1(curr_clustering_db, separator.get_cluster1()->get_dps());
			writer.write(cluster1, separator.get_cluster1()->get_assigned_scafs(), separator.get_raw_dps_cluster1(), new_cluster_name);
			cerr << '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster1.ndps() << " datapoints" << endl;
			log << '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster1.ndps() << " datapoints" << endl;
//			clustering_dbs.insert(pair<size_t, ClusterData*>(total_iterations++, cluster1));
			clusters.insert(pair<size_t, Cluster*>(total_iterations++, separator.get_cluster1()));

			sprintf(new_cluster_name, "%lu", total_iterations);
//			ClusterData* cluster2 = new ClusterData(clustering_db, separator.get_cluster2()->get_dps());
			ClusterData cluster2(curr_clustering_db, separator.get_cluster2()->get_dps());
			writer.write(cluster2, separator.get_cluster2()->get_assigned_scafs(), separator.get_raw_dps_cluster2(), new_cluster_name);
			cerr << '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster2.ndps() << " datapoints" << endl;
			log << '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster2.ndps() << " datapoints" << endl;
//			clustering_dbs.insert(pair<size_t, ClusterData*>(total_iterations++, cluster2));
			clusters.insert(pair<size_t, Cluster*>(total_iterations++, separator.get_cluster2()));
		}
		delete(&curr_clustering_db);
		delete(curr_cluster);
	}
	
	for(size_t i=0; i<=scaf_db.ndps(); i++) {
		if(dp2cluster.find(i) == dp2cluster.end())
			ofs_dp2cluster << i << "\t0" << endl; 
		else
			ofs_dp2cluster << i << "\t" << dp2cluster[i] << endl; 
	}
	for(size_t i=0; i<=scaf_db.nscafs(); i++) {
		if(scaf2cluster.find(i) == scaf2cluster.end())
			ofs_scaf2cluster << scaf_db.scaf_id2name(i) << "\t0" << endl; 
		else
			ofs_scaf2cluster << scaf_db.scaf_id2name(i) << "\t" << scaf2cluster[i] << endl; 
	}

	ofs_dp2cluster.close();
	ofs_scaf2cluster.close();

	cerr << '[' << get_time() << ']' << " Finished successfully" << endl << endl;
	log << '[' << get_time() << ']' << " Finished successfully" << endl << endl;
	log.close();
	summary.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t get_fasta(const set<size_t>& scafs, const ScafDpData& scaf_db, const char* fasta_out, const char* fasta_file)
{
	size_t total_bps = 0;

	FILE* IN = fopen(fasta_file, "r");
	FILE* OUT = fopen(fasta_out, "w");

	bool write = false;
	char line[16384];

	while(fgets(line, 16384, IN)) {
		// Option 1: this is a header line. Decide whether this entry needs to be written and write it if needed
		if(line[0] ==  '>') {
			if(write)
				fprintf(OUT, "\n");

			char *p = line;
			while(*p && !isspace(*p))
				p++;
			char c = *p;
			if(*p)
				*p = 0;
			size_t scaf_id = scaf_db.scaf_name2id(line+1);
			if(scaf_id == 0) {
				// scaffold is probably too short
				write = false;
				continue;
//				// This shouldn't happen - it probably means that the file supplied by the user is incorrect
//				cerr << __FILE__ << " (" << __LINE__ << "), fatal error: scaffold name " << (line+1) << " that appears in the assembly file was not found in the .names file" << endl << endl;
//				exit(-1);
			}
			if(scafs.find(scaf_id) == scafs.end()) {
				write = false;
			}
			else {
				write = true;
				*p = c;
				fprintf(OUT, "%s", line);

				// Just in case the header line is longer than 16,384 letters
				while((strlen(line) == 16384) && fgets(line, 16384, IN))
					fprintf(OUT, "%s", line);
			}
			continue;
		}
		// Option 2: this is part of a sequence that needs not be written or this is an empty line
		if(!write || (strlen(line) == 0))
			continue;

		// Option 3: this is part of a sequence that needs to be written
		total_bps += strlen(line);
		fprintf(OUT, "%s", line);
	}
	
	fclose(IN);
	fclose(OUT);

	return total_bps;
}

void unit_test()
{
	const char* names_file = "../testcases/Cattol.time-normalized/all.fna.no-replicates.no-zero-coverage.names";
//	const char* lrn_file = "../testcases/Cattol.time-normalized/all.fna.no-replicates.no-zero-coverage.lrn";
	const char* lrn_file = "testing/0.lrn";

	system("mkdir testing");

	cerr << '[' << get_time() << ']' << " Loading data" << endl;

	ScafDpData		scaf_db(names_file);
	ClusterData		cluster_db(lrn_file);
	ClusterEsomWriter 	writer("testing", scaf_db);
	set<size_t>		scafs;
	set<size_t>		dps;

	// 0. Write the whole file with no manipulations
	for(size_t i=1; i<=scaf_db.nscafs(); i++)
		scafs.insert(i);
	for(size_t i=1; i<=scaf_db.ndps(); i++)
		dps.insert(i);

	cerr << '[' << get_time() << ']' << " Test 0" << endl;
	writer.write(cluster_db, scafs, dps, "0");

	cerr << '[' << get_time() << ']' << " Test 1" << endl;
	// 1+2. Test the "copy constructor"
	ClusterData             cluster_db2(cluster_db, dps);
	writer.write(cluster_db2, scafs, dps, "1");

	cerr << '[' << get_time() << ']' << " Test 2" << endl;
	set<size_t>		dps3;
	for(size_t i=1; i<=scaf_db.ndps(); i+=3)
		dps3.insert(i);
	ClusterData             cluster_db3(cluster_db2, dps3);
	writer.write(cluster_db3, scafs, dps3, "2");

	// 3. Test cluster
	cerr << '[' << get_time() << ']' << " Test 3 (" << cluster_db.get_dimension(2).get_dimension_name() << ")" << endl;
	Cluster cluster1(scaf_db, cluster_db.get_dimension(2), Cluster::Le(), 0.224987);
	Cluster cluster2(scaf_db, cluster_db.get_dimension(2), Cluster::Gt(), 0.224987);

	cerr << '[' << get_time() << ']' << " Cluster 1: " << cluster1.ndps() << ", cluster2: " << cluster2.ndps() << endl;

}

