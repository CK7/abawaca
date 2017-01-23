#ifndef CLUSTERSEPARATOR_H
#define CLUSTERSEPARATOR_H

#include "ScafDpData.h"
#include "ClusterData.h"
#include "Cluster.h"
#include "Semaphore.h"
#include "SCGdb.h"
#include "ClusterQuality.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterSeparator
// A parent class for all clustering strategies.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparator {
protected:
	class ClusterStats;
	// We need to define the class here because some member functions are accessing it
	class ClusteringResult {
	public:
		ClusteringResult() : dimension(-1), value(-1)				{}
		ClusteringResult(size_t d, double v) : dimension(d), value(v)		{}
		ClusteringResult&	operator = (const ClusteringResult& cr)=delete;
		virtual void		reset()						{dimension = -1; value = -1;}
		virtual void		copy(const ClusteringResult& cr)		{dimension = cr.dimension; value = cr.value;}
		size_t			get_dimension() const				{return dimension;}
		double			get_value() const				{return value;}
		// Is this one better than cr?
		virtual bool		is_better(const ClusteringResult& cr) const = 0;
		virtual bool		is_legal() const = 0;
		virtual void		print_separation_params(ostream& os) const      {os << "Dimension: " << dimension << ", value:" << value;}
	protected:
		size_t			dimension;
		double			value;
	};
public:
	ClusterSeparator(const ScafDpData& _scaf_db, const SCGdb& _scg_db, const ClusterData& _clustering_db, Semaphore& threads_semaphore, ClusteringResult* cr) :
		scaf_db(_scaf_db), scg_db(_scg_db), cq(*(new ClusterQuality(_scaf_db, _scg_db))), clustering_db(_clustering_db), best_separation(cr), semaphore(threads_semaphore)
		{}
	// By default, separate() will iterate through all dimensions/values and call separate_specific() (see below) to do the actual decision 
	// as for which is separation is the best. It will also call initialize_specific() and post_specific() in order to do algorithm-specific
	// pre- and post- processing.
	// Returns true is the current cluster could be further split, false otherwise
	virtual bool		separate();
	Cluster*    		get_cluster1() const				{return cluster1;}
	Cluster* 	   	get_cluster2() const				{return cluster2;}
	const unordered_set<size_t>&	get_raw_dps_cluster1() const		{return raw_dps_cluster1;}
	const unordered_set<size_t>&	get_raw_dps_cluster2() const		{return raw_dps_cluster2;}
	size_t			get_separating_dimension() const		{return best_separation->get_dimension();}
	double			get_separating_value() const			{return best_separation->get_value();}
	bool			is_best_separation_legal() const		{return best_separation->is_legal();}
	virtual void		print_separation_params(ostream& os) const	{best_separation->print_separation_params(os);}
	// Move a dp from "from" to "to". Function also takes care of assigning scaffolds etc
//	virtual bool		move_dp_between_clusters(ClusterStats& from, ClusterStats& to, size_t dp) = 0;
	// Add a dp to "to", "unchanged" remains unchanged.
//	virtual bool		add_dp_to_cluster(ClusterStats& to, const ClusterStats& unchanged, size_t dp) = 0;
protected:
	// Find the best separation for one dimension
	virtual void		separate_dimension(size_t dimension_index) = 0;
	void			reconstruct_best_clusters();
	friend void 		call_separate_dimension(ClusterSeparator* cs, size_t dimension_index);
protected:
	const ScafDpData& 	scaf_db;
	const SCGdb&		scg_db;
	const ClusterData& 	clustering_db;
	ClusterQuality&		cq;
	Cluster*		cluster1 = NULL;
	Cluster*		cluster2 = NULL;
	unordered_set<size_t>	raw_dps_cluster1;
	unordered_set<size_t>	raw_dps_cluster2;
	ClusteringResult*	best_separation;
	size_t			cluster_ndps_threshold = 100;
	mutex			update_mtx;
	Semaphore&		semaphore;
};

bool comp_by_value (const Dimension::DataPoint& i, const Dimension::DataPoint& j);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparator::ClusterStats {
public:
        ClusterStats(Cluster& _cluster, const ScafDpData& _scaf_db) : scaf_db(_scaf_db), cluster(_cluster) {}
        virtual const ClusterStats&     add_dp(size_t dp) = 0;
        virtual const ClusterStats&     remove_dp(size_t dp) = 0;
	virtual const ClusterStats&	assign_scaf(size_t scaf) = 0;
	virtual const ClusterStats&	unassign_scaf(size_t scaf) = 0;
	virtual int			num_useful_ndps() const = 0;
        const Cluster&                  get_cluster() const                     {return cluster;}
protected:
        const ScafDpData&               scaf_db;
        Cluster&                        cluster;
};

#endif //CLUSTERSEPARATOR_H
