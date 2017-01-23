#include "ClusterSeparatorBySensitivitySpecificity.h"
#include <algorithm>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const ClusterSeparator::ClusterStats& ClusterSeparatorBySensitivitySpecificity::ClusterStats::add_dp(size_t dp)
{
	size_t	scaf = scaf_db.dp2scaf(dp);

	if(scaf_db.ndps(scaf) < 2)
		return *this;

	cluster += dp;

	if(cluster.assigned(scaf))
		true_positives++;
	else
		false_positives++;

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const ClusterSeparator::ClusterStats& ClusterSeparatorBySensitivitySpecificity::ClusterStats::remove_dp(size_t dp)
{
	size_t	scaf = scaf_db.dp2scaf(dp);

	if(scaf_db.ndps(scaf) < 2)
		return *this;

	cluster -= dp;

	if(cluster.assigned(scaf))
		true_positives--;
	else
		false_positives--;

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const ClusterSeparator::ClusterStats& ClusterSeparatorBySensitivitySpecificity::ClusterStats::assign_scaf(size_t scaf)
{
	size_t  scaf_total_ndps = scaf_db.ndps(scaf);

	if(cluster.assigned(scaf) || (scaf_total_ndps < 2))
		return *this;

	auto    mit = cluster.scafs_find(scaf);
	size_t  cluster_ndps = (mit == cluster.scafs_end())? 0 : mit->second;

	cluster.add_assigned_scaf(scaf);
	false_positives -= cluster_ndps;
	true_positives += cluster_ndps;
	total_dps_for_assigned_scafs += scaf_total_ndps;

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const ClusterSeparator::ClusterStats& ClusterSeparatorBySensitivitySpecificity::ClusterStats::unassign_scaf(size_t scaf)
{
	size_t  scaf_total_ndps = scaf_db.ndps(scaf);
	if(!cluster.assigned(scaf) || (scaf_total_ndps < 2))
		return *this;

	auto    mit = cluster.scafs_find(scaf);
	size_t  cluster_ndps = (mit == cluster.scafs_end())? 0 : mit->second;

	cluster.remove_assigned_scaf(scaf);
	false_positives += cluster_ndps;
	true_positives -= cluster_ndps;
	total_dps_for_assigned_scafs -= scaf_total_ndps;

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ostream& operator << (ostream& os, ClusterSeparatorBySensitivitySpecificity& cs)
{
	os << " Best separation:\t" << "dimension " << cs.get_separating_dimension() << " ("
		<< cs.clustering_db.get_dimension(cs.best_separation->get_dimension()).get_dimension_name() << ")\t"
		<< "Value " << cs.best_separation->get_value() << "\t"
		<< "Specificity " << static_cast<ClusterSeparatorBySensitivitySpecificity::ClusteringResult*>(cs.best_separation)->get_specificity() << "\t"
		<< "Sensitivity " << static_cast<ClusterSeparatorBySensitivitySpecificity::ClusteringResult*>(cs.best_separation)->get_sensitivity();

	return os;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterSeparatorBySensitivitySpecificity::separate_dimension(size_t dimension_index)
{
	semaphore.wait();
	const Dimension& dimension = clustering_db.get_dimension(dimension_index);
	const set<double>& values = dimension.get_values();
	int values_checked = 0;
	// Initialize - curr_cluster1 will have no points while curr_cluster2 will have all dps
	Cluster 						curr_cluster1(scaf_db);
	Cluster 						curr_cluster2(scaf_db);
	ClusterSeparatorBySensitivitySpecificity::ClusterStats	curr_cluster_stats1(curr_cluster1, scaf_db);
	ClusterSeparatorBySensitivitySpecificity::ClusterStats	curr_cluster_stats2(curr_cluster2, scaf_db);
	vector<Dimension::DataPoint>				dps;
	set<size_t> 						useful_scafs;

	for(vector<Dimension::DataPoint>::const_iterator vit = dimension.begin(); vit!=dimension.end(); vit++) {
		curr_cluster_stats2.add_dp(vit->dp);
		dps.push_back(*vit);
		if(scaf_db.ndps(scaf_db.dp2scaf(vit->dp)) > 1) {
			useful_scafs.insert(scaf_db.dp2scaf(vit->dp));
			curr_cluster_stats2.assign_scaf(scaf_db.dp2scaf(vit->dp));
		}
	}

	sort(dps.begin(), dps.end(), comp_by_value);
	ClusterSeparatorBySensitivitySpecificity::ClusteringResult	best_dimension_separation;

	for(auto vit=dps.begin(); vit!=dps.end();) {
		double value = vit->value;
		while(vit!=dps.end() && vit->value==value) {
			size_t scaf = scaf_db.dp2scaf(vit->dp);
			size_t scaf_total_ndps = scaf_db.ndps(scaf);

			curr_cluster_stats1.add_dp(vit->dp);
			curr_cluster_stats2.remove_dp(vit->dp);
			if(curr_cluster_stats2.get_cluster().assigned(scaf) && (scaf_total_ndps >= 2) && (2*curr_cluster_stats1.get_cluster().ndps(scaf) > scaf_total_ndps)) {
				curr_cluster_stats1.assign_scaf(scaf);
				curr_cluster_stats2.unassign_scaf(scaf);
			}
			vit++;
		}

		// Now do what we actually need to do in this specific class
		this->separate_specific(curr_cluster_stats1, curr_cluster_stats2, dimension_index, value, best_dimension_separation);
	}

	// Only a single access to the variables common to all threads
	unique_lock<mutex> lck(update_mtx);

	if(best_dimension_separation.is_better(*best_separation))
		best_separation->copy(best_dimension_separation);

	semaphore.notify();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterSeparatorBySensitivitySpecificity::separate_specific(ClusterStats& curr_cluster1, ClusterStats& curr_cluster2, size_t dimension_index, double tested_value,
	ClusterSeparatorBySensitivitySpecificity::ClusteringResult& best_dimension_separation)
{
	// Which one is smaller in terms of the number of datapoints that will be considered 
	// for the clustering?
	if((curr_cluster1.num_useful_ndps() < this->cluster_ndps_threshold) || (curr_cluster2.num_useful_ndps() < this->cluster_ndps_threshold))
		return;

	// Now pick the smaller cluster, we will compute sensitivity and specificity on that one
	const ClusterStats& small_cluster = (curr_cluster1.num_useful_ndps() < curr_cluster2.num_useful_ndps())? curr_cluster1 : curr_cluster2;
	const ClusterStats& large_cluster = (curr_cluster1.num_useful_ndps() >= curr_cluster2.num_useful_ndps())? curr_cluster1 : curr_cluster2;

	ClusterSeparatorBySensitivitySpecificity::ClusteringResult	cr(dimension_index, tested_value, small_cluster.get_sensitivity(), small_cluster.get_specificity());
	if(cr.is_better(best_dimension_separation) && cq.is_split_better(small_cluster.get_cluster(), large_cluster.get_cluster()))
		best_dimension_separation.copy(cr);
}
