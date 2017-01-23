#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include "ClusterWriter.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterWriter::write(const ClusterData& clustering_db, const set<size_t>& scafs, const unordered_set<size_t>& raw_dps, const char* cluster_name) const
{
	char	scaf_stats_file[1024] = {};
	char	scaf_file[1024] = {};

	sprintf(scaf_stats_file, "%s/%s.scaf-stats.txt", out_dir.c_str(), cluster_name);
	sprintf(scaf_file, "%s/%s.scaf-cluster.txt", out_dir.c_str(), cluster_name);

	this->write_clustering_info_file(clustering_db, cluster_name);
	this->write_scaffold_stats(clustering_db, scafs, raw_dps, scaf_stats_file);
	this->write_scaffolds(scafs, scaf_file, cluster_name);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterWriter::write_scaffolds(const set<size_t>& scafs, const char* out_file, const char* cluster_name) const
{
        FILE* fp = fopen(out_file, "w");
	if(fp == NULL) {
		cerr << "Fatal error, " << __FILE__ << " (" << __LINE__ << "): could write to file " << out_file << endl << endl;
		exit(-1);
	}

	for(auto vit = scafs.begin(); vit != scafs.end(); vit++)
		fprintf(fp, "%s\t%s\n", scaf_db.scaf_id2name(*vit).c_str(), cluster_name);
	fclose(fp);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterWriter::write_scaffold_stats(const ClusterData& clustering_db, const set<size_t>& scafs, const unordered_set<size_t>& raw_dps, const char* out_file) const
{
	map<size_t, size_t>	scafs_in2ndps, scafs_out2ndps;

	for(auto vit=scafs.begin(); vit != scafs.end(); vit++)
		scafs_in2ndps.insert(pair<size_t, size_t>(*vit, 0));

	for(auto vit=raw_dps.begin(); vit != raw_dps.end(); vit++) {
		size_t scaf = scaf_db.dp2scaf(*vit);
		if(scafs_in2ndps.find(scaf) != scafs_in2ndps.end())
			scafs_in2ndps[scaf]++;
		else if(scafs_out2ndps.find(scaf) != scafs_out2ndps.end())
			scafs_out2ndps[scaf]++;
		else
			scafs_out2ndps.insert(pair<size_t, size_t>(scaf, 1));
	}

        FILE* fp = fopen(out_file, "w");
	if(fp == NULL) {
		cerr << "Fatal error, " << __FILE__ << " (" << __LINE__ << "): could write to file " << out_file << endl << endl;
		exit(-1);
	}

	fprintf(fp, "SCAF      \tRAW\tTOTAL\tIN/OUT_CLUSTER\n\n");

	for(map<size_t, size_t>::const_iterator mit = scafs_in2ndps.begin(); mit != scafs_in2ndps.end(); mit++)
		fprintf(fp, "%s\t%lu\t%lu\tin\n", scaf_db.scaf_id2name(mit->first).c_str(), mit->second, scaf_db.ndps(mit->first));
	for(map<size_t, size_t>::const_iterator mit = scafs_out2ndps.begin(); mit != scafs_out2ndps.end(); mit++)
		fprintf(fp, "%s\t%lu\t%lu\tout\n", scaf_db.scaf_id2name(mit->first).c_str(), mit->second, scaf_db.ndps(mit->first));
	fclose(fp);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterEsomWriter::write_clustering_info_file(const ClusterData& clustering_db, const char* cluster_name) const
{
	char	lrn_file[1024] = {};

	sprintf(lrn_file, "%s/%s.lrn", out_dir.c_str(), cluster_name);

        FILE* fp = fopen(lrn_file, "w");
	if(fp == NULL) {
		cerr << "Fatal error, " << __FILE__ << " (" << __LINE__ << "): could write to file " << lrn_file << endl << endl;
		exit(-1); 
	}

        fprintf(fp, "%c %lu\n", '%', clustering_db.ndps());
        fprintf(fp, "%c %lu\n", '%', clustering_db.ndimensions()+1);
        fprintf(fp, "%c 9", '%');
        for(size_t i=1; i<=clustering_db.ndimensions(); i++)
                fprintf(fp, "\t1");
        fprintf(fp, "\n");
        fprintf(fp, "%c Key", '%');
        for(size_t i=1; i<=clustering_db.ndimensions(); i++)
                fprintf(fp, "\t%s", clustering_db.get_dimension(i).get_dimension_name().c_str());
        fprintf(fp, "\n");

	const unordered_set<size_t>& dps = clustering_db.datapoints();
        for(auto vit=dps.begin(); vit!=dps.end(); vit++) {
                fprintf(fp, "%lu", *vit);
                for(size_t i=1; i<=clustering_db.ndimensions(); i++)
                        fprintf(fp, "\t%lf", clustering_db.get_dimension(i).get_value(*vit));
                fprintf(fp, "\n");
        }
        fclose(fp);
}

