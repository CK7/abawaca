#ifndef CLUSTERSEPARATOR_H
#define CLUSTERSEPARATOR_H

#include "ScafDpData.h"
#include "ClusterData.h"
#include "Cluster.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterSeparator
// A parent class for all clustering strategies.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparator {
public:
	ClusterSeparator(const ScafDpData& _scaf_db, const ClusterData& _clustering_db) : scaf_db(_scaf_db), clustering_db(_clustering_db),
		cluster1(NULL), cluster2(NULL), separating_dimension(-1), separating_value(-1)	
		{}
	// By default, separate() will iterate through all dimensions/values and call separate_specific() (see below) to do the actual decision 
	// as for which is separation is the best. It will also call initialize_specific() and post_specific() in order to do algorithm-specific
	// pre- and post- processing.
	// Returns true is the current cluster could be further split, false otherwise
	virtual bool		separate();
	Cluster*    		get_cluster1() const			{return cluster1;}
	Cluster* 	   	get_cluster2() const			{return cluster2;}	
	const set<size_t>&	get_raw_dps_cluster1() const		{return raw_dps_cluster1;}
	const set<size_t>&	get_raw_dps_cluster2() const		{return raw_dps_cluster2;}
	size_t			get_separating_dimension() const	{return separating_dimension;}
	double			get_separating_value() const		{return separating_value;}
protected:
	// This function will be called after initialization of ClusterSeparator's protected fields whem ClusterSeparator::separate starts to run
	virtual void		initialize_specific()			{}
	// This function will be called while ClusterSeparator::separate iterates through the different dimensions/values in order to check each separation
	virtual void 		separate_specific(Cluster& curr_cluster1, Cluster& curr_cluster2, size_t cluster1_useful_ndps, size_t cluster2_useful_ndps, size_t dimension_index, double tested_value) = 0;
	// This function will be the very last piece of code ClusterSeparator::separate executes. its return value is also ClusterSeparator::separate return value.
	virtual bool		post_specific()				{return true;}
protected:
	const ScafDpData& 	scaf_db;
	const ClusterData& 	clustering_db;
	Cluster*		cluster1;
	Cluster*		cluster2;
	set<size_t>		raw_dps_cluster1;
	set<size_t>		raw_dps_cluster2;
	size_t			separating_dimension;
	double			separating_value;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterSeparatorBySensitivitySpecificity
// Implements the basic strategy that was successfully attempted in the initial PERL version of this algorithm
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparatorBySensitivitySpecificity : public ClusterSeparator {
public:
	ClusterSeparatorBySensitivitySpecificity(const ScafDpData& scaf_db, const ClusterData& clustering_db) : ClusterSeparator(scaf_db, clustering_db),
		sensitivity_threshold(0.8), specificity_threshold(0.8), prodduct_threshold(0.80), sum_threshold(1.6), best_sensitivity(0), 
		best_specificity(0)
		{}
	double			get_best_sensitivity() const		{return best_sensitivity;}
	double			get_best_specificity() const		{return best_specificity;}
protected:
	void 			separate_specific(Cluster& curr_cluster1, Cluster& curr_cluster2, size_t cluster1_useful_ndps, size_t cluster2_useful_ndps, size_t dimension_index, double tested_value);
	bool			legal(double sensitivity, double specificity) const;
	bool			best(double sensitivity, double specificity) const;
protected:
	double			sensitivity_threshold;
	double			specificity_threshold;
	double			prodduct_threshold;
	double			sum_threshold;
	double			best_sensitivity, best_specificity;
};

#endif //CLUSTERSEPARATOR_H
