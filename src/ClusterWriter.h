#ifndef CLUSTERWRITER_H
#define CLUSTERWRITER_H

#include "ScafDpData.h"
#include "ClusterData.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterWriter
// This class serves as a parent class for all classes that have to do with writing cluster information. At the moment there is just 
// one decsendent class - ClusterEsomWriter - that outputs file in the Databionic ESOM .lrn file. All descendent classes have to 
// output the following 3 types of information:
// 1. clustering data - for each datapoint in the cluster, this file should contain the data that was used for the clustering.
// 2. Scaffolds file - a list of all scaffolds that belong to the current cluster
// 3. scaffold statistics - number of datapoints from each scaffold that were included in the original clustering, before any 
//    post-clustering manipulations (e.g. adding datapoints that were not originally in the cluster based on their relation to 
//    scaffolds that belong to the cluster)   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterWriter {
public:
	ClusterWriter(string _out_dir, const ScafDpData& _scaf_db) : out_dir(_out_dir), scaf_db(_scaf_db) 		{}
	void 		write(const ClusterData& clustering_db, const set<size_t>& scafs, const unordered_set<size_t>& raw_dps, const char* cluster_name) const;
protected:
	void		write_scaffolds(const set<size_t>& scafs, const char* out_file, const char* cluster_name) const;
	void		write_scaffold_stats(const ClusterData& clustering_db, const set<size_t>& scafs, const unordered_set<size_t>& raw_dps, const char* out_file) const;
	virtual void	write_clustering_info_file(const ClusterData& clustering_db, const char* cluster_name) const = 0;
protected:
	string	 		out_dir;
	const ScafDpData& 	scaf_db;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterEsomWriter : public ClusterWriter {
public:
	ClusterEsomWriter(string _out_dir, const ScafDpData& _scaf_db) : ClusterWriter(_out_dir, _scaf_db)		{}
//	void 		write(const ClusterData& clustering_db, const set<size_t>& scafs, const set<size_t>& raw_dps, const char* cluster_name) const;
protected:
	void	write_clustering_info_file(const ClusterData& clustering_db, const char* cluster_name) const;
};

#endif // CLUSTERWRITER_H
