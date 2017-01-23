#ifndef CLUSTERQUALITY_H
#define CLUSTERQUALITY_H

#include "Cluster.h"
#include "ScafDpData.h"
#include "SCGdb.h"
#include "Dimension.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Visitor class to Cluster, performs operations that can be used to determine the quality of the cluster
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterQuality {
public:
	ClusterQuality(const ScafDpData& _scaf_db, const SCGdb& _scg_db) : scaf_db(_scaf_db), scg_db(_scg_db)	{}
	virtual 		~ClusterQuality()									{}
	void			gc(const Cluster& cluster, double& avg_gc, double& gc_stdev) const;
	void			cvg(const Cluster& cluster, double& avg_cvg, double& cvg_stdev) const;
	bool			is_split_better(/*const Cluster& cluster, */const Cluster& sub_cluster1, const Cluster& sub_cluster2) const;
	void			scg(const Cluster& cluster, double& nunique, double& avg) const;
	void			dimension(const Cluster& cluster, double& avg, double& stdev, const Dimension& dimension) const;
	size_t 			num_scaffolds(const Cluster& cluster) const;
	size_t			total_size(const Cluster& cluster) const;
protected:
	const ScafDpData&	scaf_db;
	const SCGdb&		scg_db;
};

#endif //CLUSTERQUALITY_H
