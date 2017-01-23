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
#include "ClusterSeparatorBySensitivitySpecificity.h"
#include "ClusterSeparatorSplitScafs.h"
#include "ClusterWriter.h"
#include "ClusterQuality.h"
#include "Semaphore.h"
#include "SCGdb.h"
#include "SeqIOWrite_fasta.h"

using namespace std;

#define	VERSION "v1.07"

size_t get_fasta(const set<size_t>& scafs, const ScafDpData& scaf_db, const char* fasta_out);
Cluster* init_cluster(const ScafDpData& scaf_db);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Prog_params {
	int 			read_params(int argc, const char** argv);
	void 			usage(const char* prog_name) const;
	static Prog_params*	Instance();
protected:
	Prog_params()	{}
	void	read_from_build_dir(int argc, const char** argv);
public:
	string	lrn_file;
	string	names_file;
	string	info_file;
	string	links_file;
	string	fasta_file;
	string	root_dir;
	string	gene2scg_file;
	string	scg_list_file = "/home/itaish/software/cpp/abawaca/curr-version/scg.list";
	string	abawaca_build_dir;
	string	cluster_dir;
	string	fasta_dir;
	int	ncpus = 1;
public:
	static Prog_params*	_instance;
};

Prog_params* Prog_params::_instance = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char** argv) 
{
	cerr << "abawaca " << VERSION << endl << endl;

	Prog_params&		params = *Prog_params::Instance();

	ofstream		log;
	ofstream		summary;
	ofstream		ofs_dp2cluster;
	ofstream		ofs_scaf2cluster;
	char 			line[1024] = {};

	if(params.read_params(argc, argv))
		return -1;

	Semaphore		thread_semphore(params.ncpus);

	sprintf(line, "%s/log", params.root_dir.c_str());
	log.open(line);
	sprintf(line, "%s/summary.txt", params.root_dir.c_str());
	summary.open(line);
	sprintf(line, "%s/dp2cluster.txt", params.root_dir.c_str());
	ofs_dp2cluster.open(line);
	sprintf(line, "%s/scaf2cluster.txt", params.root_dir.c_str());
	ofs_scaf2cluster.open(line);

	cerr << '[' << get_time() << ']' << " Starting" << endl;
	cerr << '[' << get_time() << ']' << " Creating databases based on " << params.names_file << " and " << params.lrn_file << endl;
	log << "abawaca " << VERSION << endl << endl;
	log << '[' << get_time() << ']' << " Starting" << endl;
	log << '[' << get_time() << ']' << " Creatinf databases based on " << params.names_file << " and " << params.lrn_file << endl;
	ScafDpData			scaf_db(params.names_file, params.fasta_file, params.info_file);
	ClusterData			clustering_db(params.lrn_file, scaf_db);
	SCGdb				scg_db(scaf_db, params.gene2scg_file, params.scg_list_file);
	map<size_t, Cluster*>		clusters;
	size_t				total_iterations = 1;
	map<size_t, size_t>		dp2cluster;
	map<size_t, size_t>		scaf2cluster;

	summary << "Cluster\t# scafs\t# dps\t# bps\t%G+C\tStdev\tCvg\tstdev\t#SCG" << '/' << scg_db.get_total_num_scgs() << "\tAvg" << endl;

	clusters[total_iterations++] = init_cluster(scaf_db);

	// We will keep iterating over all items in clustering_dbs until it is empty.
	while(clusters.size() > 0) {
		map<size_t, Cluster*>::iterator 	mit = clusters.begin();
		Cluster* 				curr_cluster = mit->second;
		ClusterData&    			curr_clustering_db = *(new ClusterData(clustering_db, curr_cluster->get_dps()));
		size_t 					iteration = mit->first;

		clusters.erase(mit);

		cerr << endl << '[' << get_time() << ']' << " Checking cluster " << iteration << " (" << curr_cluster->ndps() << " datapoints)" << endl;
		log << endl << '[' << get_time() << ']' << " Checking cluster " << iteration << " (" << curr_cluster->ndps() << " datapoints)" << endl;

		ClusterSeparatorBySensitivitySpecificity separator(scaf_db, scg_db, curr_clustering_db, thread_semphore);
//		ClusterSeparatorSplitScafs separator(scaf_db, scg_db, curr_clustering_db, thread_semphore);

		separator.separate();

		cerr << '[' << get_time() << ']' << " ";
		separator.print_separation_params(cerr);
		cerr << endl;
		log << '[' << get_time() << ']' << " ";
		separator.print_separation_params(log);
		log << endl;

		ClusterQuality cq(scaf_db, scg_db);
		double	nunique_scgs, avg_ncopies_scgs;
		double	avg_gc, stdev_gc;
		double	avg_cvg=-1, stdev_cvg=-1;

		if((separator.get_cluster1() == NULL)/* || (!cq.is_split_better(*(separator.get_cluster1()), *(separator.get_cluster2())))*/) {
			if(separator.get_cluster1() == NULL) {
				cerr << '[' << get_time() << ']' << " Cluster " << iteration << " is a terminal cluster" << endl;
				log << '[' << get_time() << ']' << " Cluster " << iteration << " is a terminal cluster" << endl;
			}
			else {
				cerr << '[' << get_time() << ']' << " Cluster " << iteration << " will not be split further (terminal cluster)" << endl;
				log << '[' << get_time() << ']' << " Cluster " << iteration << " will not be split further (terminal cluster)" << endl;
			}
			for(auto sit = curr_cluster->dps_begin(); sit != curr_cluster->dps_end(); sit++)
				dp2cluster[*sit] = iteration;
			for(auto sit = curr_cluster->assigned_scafs_begin(); sit != curr_cluster->assigned_scafs_end(); sit++)
				scaf2cluster[*sit] = iteration;
			char fasta_out[1024];

			sprintf(fasta_out, "%s/%lu.fasta", params.fasta_dir.c_str(), iteration);
			size_t total_bps = get_fasta(curr_cluster->get_assigned_scafs(), scaf_db, fasta_out);

			total_bps = cq.total_size(*curr_cluster);
			size_t	nscaffolds = cq.num_scaffolds(*curr_cluster);
			cq.scg(*curr_cluster, nunique_scgs, avg_ncopies_scgs);
			cq.gc(*curr_cluster,  avg_gc, stdev_gc);
			cq.cvg(*curr_cluster,  avg_cvg, stdev_cvg);

			summary << iteration << '\t' << nscaffolds << '\t' << curr_cluster->ndps() << '\t' << total_bps << '\t'
				<< int(1000*avg_gc)/10.0 << '\t' << int(1000*stdev_gc)/100.0 << '\t' << int(10*avg_cvg)/10.0
				<< '\t' << int(10*stdev_cvg)/10.0 << '\t' << nunique_scgs << '\t' << int(100*avg_ncopies_scgs)/100.0 << endl;

		}
		else {
			cerr << '[' << get_time() << ']' << separator << endl;
			cerr << '[' << get_time() << ']' << " Cluster " << iteration << " will be separated to the following two clusters" << endl;
			log << '[' << get_time() << ']' << separator << endl;
			log << '[' << get_time() << ']' << " Cluster " << iteration << " will be separated to the following two clusters" << endl;

			ClusterEsomWriter writer(params.cluster_dir, scaf_db);
			char new_cluster_name[128];

			sprintf(new_cluster_name, "%lu", total_iterations);
			ClusterData cluster1(curr_clustering_db, separator.get_cluster1()->get_dps());
			writer.write(cluster1, separator.get_cluster1()->get_assigned_scafs(), separator.get_raw_dps_cluster1(), new_cluster_name);

			cq.scg(*(separator.get_cluster1()), nunique_scgs, avg_ncopies_scgs);
			cq.gc(*(separator.get_cluster1()), avg_gc, stdev_gc);
			cq.cvg(*(separator.get_cluster1()), avg_cvg, stdev_cvg);

			cerr 	<< '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster1.ndps() << " datapoints, "
				<< nunique_scgs << " unique SCGs (" << avg_ncopies_scgs << " copies), %G+C=" << avg_gc << ", coverage=" << avg_cvg << endl;
			log 	<< '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster1.ndps() << " datapoints, "
				<< nunique_scgs << " unique SCGs (" << avg_ncopies_scgs << " copies), %G+C=" << avg_gc << ", coverage=" << avg_cvg << endl;
//			clustering_dbs.insert(pair<size_t, ClusterData*>(total_iterations++, cluster1));
			clusters.insert(pair<size_t, Cluster*>(total_iterations++, separator.get_cluster1()));

			sprintf(new_cluster_name, "%lu", total_iterations);
//			ClusterData* cluster2 = new ClusterData(clustering_db, separator.get_cluster2()->get_dps());
			ClusterData cluster2(curr_clustering_db, separator.get_cluster2()->get_dps());
			writer.write(cluster2, separator.get_cluster2()->get_assigned_scafs(), separator.get_raw_dps_cluster2(), new_cluster_name);

			cq.scg(*(separator.get_cluster2()), nunique_scgs, avg_ncopies_scgs);
			cq.gc(*(separator.get_cluster2()), avg_gc, stdev_gc);
			cq.cvg(*(separator.get_cluster2()), avg_cvg, stdev_cvg);

			cerr 	<< '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster2.ndps() << " datapoints, "
				<< nunique_scgs << " unique SCGs (" << avg_ncopies_scgs << " copies), %G+C=" << avg_gc << ", coverage=" << avg_cvg << endl;
			log 	<< '[' << get_time() << ']' << " Cluster " << total_iterations << ", " << cluster2.ndps() << " datapoints, "
				<< nunique_scgs << " unique SCGs (" << avg_ncopies_scgs << " copies), %G+C=" << avg_gc << ", coverage=" << avg_cvg << endl;
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
	for(size_t i=1; i<=scaf_db.nscafs(); i++) {
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
size_t get_fasta(const set<size_t>& scafs, const ScafDpData& scaf_db, const char* fasta_out)
{
	size_t total_bps = 0;
	Bio::SeqIOWrite_fasta<Bio::DNASequence>  writer(fasta_out);


	for(auto sit=scafs.begin(); sit!=scafs.end(); sit++) {
		const Bio::DNASequence* seq_obj = scaf_db.get_scaf(*sit)->get_seq_obj();
		if(seq_obj == NULL) {
			cerr << "Warning: could not find scaffold " << scaf_db.scaf_id2name(*sit) << " (" << *sit << ") in fasta file" << endl;
			continue;
		}
		writer.write_seq(*seq_obj);
		total_bps += seq_obj->seq().size();
	}

	return total_bps;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cluster* init_cluster(const ScafDpData& scaf_db)
{
	Cluster* cluster = new Cluster(scaf_db);
	for(size_t i=1; i<=scaf_db.ndps(); i++)
		*cluster += i;
	for(size_t i=1; i<=scaf_db.nscafs(); i++)
		cluster->add_assigned_scaf(i);

	return cluster;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Prog_params* Prog_params::Instance()
{
	if(Prog_params::_instance == NULL)
		Prog_params::_instance = new Prog_params();

	return Prog_params::_instance;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Prog_params::usage(const char* prog_name) const
{
	fprintf(stderr, "Usage: %s -u <abawaca-build-directory> -f <fasta-file> -o <out-dir> [-p <# CPUs>] [-c <gene2scg>]\n\n", prog_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Prog_params::read_params(int argc, const char** argv)
{
	if((argc == 1) || ((argc == 2) && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")))) {
		usage(argv[0]);
		return -1;
	}

	// Read parameters
	for(auto i=1; i<argc; i++) {
		if(!strcmp(argv[i], "-u")) {
			abawaca_build_dir = argv[++i];
			char data_file[2048];
			sprintf(data_file, "%s/data.txt", abawaca_build_dir.c_str());
			FILE* fp = fopen(data_file, "r");
			if(fp == NULL) {
				cerr << "Error: could not find " << abawaca_build_dir << "/data.txt. Is it possible that the data was generated using an old version of abawaca-build? (prior to v1.07)" << endl << endl; 
				return -1;
			}
			char line[8192];
			while(fgets(line, 8192, fp)) {
				char* p = strrchr(line, '\n');
				while((p != NULL) && (*p == '\n'))  {
					*(p--) = 0;
				}
				p=strchr(line, '\t');
				*(p++) = 0;
				while(isspace(*p) && (*p != 0))
					p++;

				if(!strcmp(line, "Links")) {
					links_file = p;
				}
				else if(!strcmp(line, "Info")) {
					info_file = p;
				}
				else if(!strcmp(line, "Names")) {
					names_file = p;
				}
				else if(!strcmp(line, "Lrn")) {
					lrn_file = p;
				}
				else if(!strcmp(line, "Assembly")) {
					fasta_file = p;
				}
				else if(!strcmp(line, "SCG")) {
					gene2scg_file = p;
				}
			}
/*			lrn_file = abawaca_build_dir + "/abawaca.lrn";
			names_file = abawaca_build_dir + "/abawaca.names";
			info_file = abawaca_build_dir + "/abawaca.info";
			links_file = abawaca_build_dir + "/abawaca.links";*/
		}
		else if(!strcmp(argv[i], "-f"))
			fasta_file = argv[++i];
		else if(!strcmp(argv[i], "-o"))
			root_dir = argv[++i];
		else if(!strcmp(argv[i], "-p")) {
			ncpus = atoi(argv[++i]);
			if(ncpus==0 || ncpus>40) {
				cerr << "Error: # CPUs (" << argv[i] << ") must be an integer between 1 to 40" << endl << endl;
				return -1; 
			}
		}
		else if(!strcmp(argv[i], "-c"))
			gene2scg_file = argv[++i];
		else {
			cerr << "Unrecognizes flag: " << argv[i] << endl;
			return -1;
		}
	}

	// Make sure we have everything we need
	if(lrn_file == "") {
		cerr << "abawaca-build directrory was not specified (-u)" << endl << endl;
		return -1;
	}

	FILE* IN = fopen(names_file.c_str(), "r");
	if(IN == NULL) {
		cerr << "Could not read names file (" << names_file << ") from the abawaca directory (" << abawaca_build_dir << ")" << endl;
		return -1;
	}
	fclose(IN);

	IN = fopen(lrn_file.c_str(), "r");
	if(IN == NULL) {
		cerr << "Could not read lrn file (" << lrn_file << ") from the abawaca directory (" << abawaca_build_dir << ")" << endl;
		return -1;
	}
	fclose(IN);

	IN = fopen(info_file.c_str(), "r");
	if(IN == NULL) {
		cerr << "Could not read info file (" << info_file << ") from the abawaca directory (" << abawaca_build_dir << ")" << endl;
		return -1;
	}
	fclose(IN);

	IN = fopen(links_file.c_str(), "r");
	if(IN == NULL) {
		cerr << "Could not read links file (" << links_file << ") from the abawaca directory (" << abawaca_build_dir << ")" << endl;
		return -1;
	}
	fclose(IN);

	IN = fopen(fasta_file.c_str(), "r");
	if(IN == NULL) {
		cerr << "Could not read " << fasta_file << endl;
		return -1;
	}
	fclose(IN);

	if(gene2scg_file.size() > 0) {
		IN = fopen(gene2scg_file.c_str(), "r");
		if(IN == NULL) {
			cerr << "Could not read " << gene2scg_file << endl;
			return -1;
		}
		fclose(IN);
	}

//	const char* p = strrchr(argv[0], '/');
//	if(p != NULL) {
//		string s = argv[0];
//		scg_list_file = s.substr(0, p-argv[0]+1) + scg_list_file;
//	}
	IN = fopen(scg_list_file.c_str(), "r");
	if(IN == NULL) {
		cerr << "Could not read " << scg_list_file << endl;
		return -1;
	}
	fclose(IN);

	cluster_dir = root_dir + string("/clusters");
	fasta_dir = root_dir + string("/final-clusters");

	char 			line[1024] = {};

	if(!directory_exists(root_dir.c_str())) {
		sprintf(line, "mkdir %s", root_dir.c_str());
		system(line);
	}
	if(!directory_exists(cluster_dir.c_str())) {
		sprintf(line, "mkdir %s", cluster_dir.c_str());
		system(line);
	}
	if(!directory_exists(fasta_dir.c_str())) {
		sprintf(line, "mkdir %s", fasta_dir.c_str());
		system(line);
	}

	return 0;
}
