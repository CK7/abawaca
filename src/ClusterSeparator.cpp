#include "ClusterSeparator.h"
#include <iostream>
#include <thread>
#include <exception>
#include "misc.h"
#include "Semaphore.h"

bool comp_by_value (const Dimension::DataPoint& i, const Dimension::DataPoint& j) { return (i.value<j.value); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClusterSeparator::ClusteringResult::is_better(const ClusteringResult& cr) const
{
        if(dimension != cr.dimension) 
                return dimension < cr.dimension; 
        return value < cr.value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void call_separate_dimension(ClusterSeparator* cs, size_t dimension_index)
{
	cs->separate_dimension(dimension_index);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterSeparator::reconstruct_best_clusters()
{
	if(!this->is_best_separation_legal())
		return;

	const Dimension& 		dimension = clustering_db.get_dimension(this->get_separating_dimension());

	this->cluster1 = new Cluster(scaf_db);
	this->cluster2 = new Cluster(scaf_db);

	size_t cluster1_useful_ndps=0, cluster2_useful_ndps=0;
	for(vector<Dimension::DataPoint>::const_iterator vit = dimension.begin(); vit!=dimension.end(); vit++) {
		int useful = (scaf_db.ndps(scaf_db.dp2scaf(vit->dp)) > 1)? 1 : 0;

		if(vit->value <= this->get_separating_value()) {
			*(this->cluster1) += vit->dp;
			cluster1_useful_ndps += useful;
		}
		else {
			*(this->cluster2) += vit->dp;
			cluster2_useful_ndps += useful;
		}
	}

	if(cluster2_useful_ndps < cluster1_useful_ndps) {
		Cluster* temp = this->cluster1;
		this->cluster1 = this->cluster2;
		this->cluster2 = temp;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClusterSeparator::separate()
{
	// Start by cleaning everything from previous runs
	if(this->cluster1 != NULL)
		delete(this->cluster1);
	if(this->cluster2 != NULL)
		delete(this->cluster2);
	this->cluster1 = this->cluster2 = NULL;
	this->best_separation->reset();
	this->raw_dps_cluster1.clear();
	this->raw_dps_cluster2.clear();

/*	for(size_t dimension_index=1; dimension_index<=clustering_db.ndimensions(); dimension_index++) {
		this->separate_dimension(dimension_index);
	}*/
	thread*	threads = new thread[clustering_db.ndimensions()];
	for(size_t dimension_index=1; dimension_index<=clustering_db.ndimensions(); dimension_index++) {
		threads[dimension_index-1] = thread(call_separate_dimension, this, dimension_index);
	}
	for(auto i=0; i<clustering_db.ndimensions(); i++)
		threads[i].join();
	delete[] (threads);

	this->reconstruct_best_clusters();

	// If no significant separation was found simply return false
	if(cluster1 == NULL)
		return false;

	// Before switching datapoints between the clusters keep the original sets
	for(auto sit = cluster1->dps_begin(); sit != cluster1->dps_end(); sit++)
		raw_dps_cluster1.insert(*sit);
	for(auto sit = cluster2->dps_begin(); sit != cluster2->dps_end(); sit++)
		raw_dps_cluster2.insert(*sit);

	// Now decide, for each scaffold, to which cluster it belongs. Cluster1 will contain all scaffolds that have at least half
	// of their datapoints in it, cluster2 will contain the rest
	for(auto mit = cluster1->scafs_begin(); mit != cluster1->scafs_end(); mit++)
		if(2*mit->second >= scaf_db.ndps(mit->first)) {
			cluster1->add_assigned_scaf(mit->first);
		}
	for(auto mit = cluster2->scafs_begin(); mit != cluster2->scafs_end(); mit++)
		if(2*mit->second > scaf_db.ndps(mit->first)) {
			cluster2->add_assigned_scaf(mit->first);
		}

	// Now move the datapoints based on their scaffold assignment
	set<size_t> dps;
	for(auto sit = cluster1->dps_begin(); sit != cluster1->dps_end(); sit++)
		if(!cluster1->assigned(scaf_db.dp2scaf(*sit)))
			dps.insert(*sit);

	for(auto sit = dps.begin(); sit != dps.end(); sit++) {
		*cluster1 -= *sit;
		*cluster2 += *sit;
	}

	dps.clear();
	for(auto sit = cluster2->dps_begin(); sit != cluster2->dps_end(); sit++)
		if(!cluster2->assigned(scaf_db.dp2scaf(*sit)))
			dps.insert(*sit);

	for(auto sit = dps.begin(); sit != dps.end(); sit++) {
		*cluster1 += *sit;
		*cluster2 -= *sit;
	}

	// Check again if the resulting clusters are sufficiently large.
	if((cluster1->ndps() < this->cluster_ndps_threshold) || (cluster2->ndps() < this->cluster_ndps_threshold)) {
		delete(cluster1);
		cluster1 = NULL;
		delete(cluster2);
		cluster2 = NULL;
		best_separation->reset();
		return false;
	}

	return true;
}
