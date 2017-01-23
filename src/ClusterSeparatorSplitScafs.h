#ifndef CLUSTERSEPARATORSPLITSCAFS_H
#define CLUSTERSEPARATORSPLITSCAFS_H

#include <exception>
#include "ClusterSeparator.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterSeparatorSplitScafs
// Implements the following strategy:
// 1. Find how many scaffolds are split. A scaffold is split if both clusters contain less than fraction_dps_to_be_considered_in of its dps.
// 2. Calculate the number of scaffolds in the smaller cluster
// 3. Look for the ratio between the number of split scaffolds and the smaller cluster. if this ratio is >= split_scaf_ratio_threshold
//    then the split is legal and can be used. 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparatorSplitScafs : public ClusterSeparator {
protected:
	class ClusterStats;
	struct IsScafIn;
	class ClusteringResult : public ClusterSeparator::ClusteringResult{
	public:
		ClusteringResult() : ClusterSeparator::ClusteringResult()	{}
		ClusteringResult(size_t d, double v, double s, double c) : ClusterSeparator::ClusteringResult(d, v), split_scaf_ratio(s), cluster_size_ratio(c)	{}
		void			reset() override				{split_scaf_ratio = 1000000; cluster_size_ratio = 1000000; ClusterSeparator::ClusteringResult::reset();}
		void			copy(const ClusterSeparator::ClusteringResult& cr) override
		{
			this->split_scaf_ratio = dynamic_cast<const ClusteringResult&>(cr).split_scaf_ratio;
			this->cluster_size_ratio = dynamic_cast<const ClusteringResult&>(cr).cluster_size_ratio;
			ClusterSeparator::ClusteringResult::copy(cr);
		}
		// Is this one better than cr?
		virtual bool		is_better(const ClusterSeparator::ClusteringResult& cr_) const override
		{
			const ClusteringResult& cr = dynamic_cast<const ClusteringResult&>(cr_);
			if(split_scaf_ratio != cr.split_scaf_ratio) return split_scaf_ratio < cr.split_scaf_ratio;
			if(cluster_size_ratio != cr.cluster_size_ratio) return cluster_size_ratio < cr.cluster_size_ratio;
			return ClusterSeparator::ClusteringResult::is_better(cr);
		}
		bool			is_legal() const override			{return (split_scaf_ratio <= split_scaf_ratio_threshold);}
		double			get_split_scaf_ratio() const			{return split_scaf_ratio;}
		double			get_cluster_size_ratio() const			{return cluster_size_ratio;}
	protected:
		double			split_scaf_ratio = 1000000;
		double			cluster_size_ratio = 1000000;
		const double		split_scaf_ratio_threshold = 0.1;
	};
public:
	ClusterSeparatorSplitScafs(const ScafDpData& scaf_db, const SCGdb& scg_db, const ClusterData& clustering_db,  Semaphore& threads_semaphore) : 
		ClusterSeparator(scaf_db, scg_db, clustering_db, threads_semaphore, new ClusteringResult())
		{}
	double			get_best_split_scaf_ratio() const	{return dynamic_cast<const ClusteringResult&>(*best_separation).get_split_scaf_ratio();}
	friend ostream& 	operator << (ostream& os, ClusterSeparatorSplitScafs& cs);
protected:
	void			separate_dimension(size_t dimension_index);
	void 			separate_specific(ClusterStats& curr_cluster1, ClusterStats& curr_cluster2, size_t dimension_index, double tested_value, size_t nuseful_scafs,
						  ClusterSeparatorSplitScafs::ClusteringResult& best_dimension_separation);
	void 			separate_specific(Cluster& curr_cluster1, Cluster& curr_cluster2, size_t cluster1_useful_ndps, size_t cluster2_useful_ndps, size_t dimension_index, double tested_value)
		{throw std::exception();}
protected:
	double			fraction_dps_to_be_considered_in = 0.8;
	double			split_scaf_ratio_threshold = 0.1;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparatorSplitScafs::ClusterStats {
public:
	ClusterStats(Cluster& _cluster, const ScafDpData& _scaf_db, IsScafIn& sb) : scaf_db(_scaf_db), cluster(_cluster), is_scaf_part_of(sb) {}
	const ClusterStats&	operator += (size_t dp);
	const ClusterStats&	operator -= (size_t dp);
	int			ndps_this_cluster_scafs() const		{return ndps_in_scafs_that_belong;}
	int			num_useful_ndps() const			{return cluster_useful_ndps;}
	int			num_scafs_in_this_cluster() const	{return scafs_in;}
	const Cluster&		get_cluster() const			{return cluster;}
protected:
	const ScafDpData&	scaf_db;
	Cluster& 		cluster;
	int			cluster_useful_ndps = 0;
	int			ndps_in_scafs_that_belong = 0;
	int			scafs_in = 0;
	IsScafIn&		is_scaf_part_of;
};

struct ClusterSeparatorSplitScafs::IsScafIn {
	IsScafIn(double _fraction_dps_to_be_considered_in) : fraction_dps_to_be_considered_in(_fraction_dps_to_be_considered_in)
		{}
	bool operator () (size_t dps_in_cluster, size_t dps_in_scaf) const
		{return (dps_in_cluster >= fraction_dps_to_be_considered_in*dps_in_scaf);}
protected:
	double			fraction_dps_to_be_considered_in;
};

ostream& operator << (ostream& os, ClusterSeparatorSplitScafs& cs);

#endif //CLUSTERSEPARATORSPLITSCAFS_H
