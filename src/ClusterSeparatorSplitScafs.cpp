#include "ClusterSeparatorSplitScafs.h"
#include <iostream>
#include <algorithm>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const ClusterSeparatorSplitScafs::ClusterStats& ClusterSeparatorSplitScafs::ClusterStats::operator += (size_t dp)
{
	auto	mit = cluster.scafs_find(this->scaf_db.dp2scaf(dp));
	size_t	ndps_in_cluster = (mit != cluster.scafs_end())? mit->second : 0;
	size_t	ndps_in_scaf = scaf_db.ndps(this->scaf_db.dp2scaf(dp));

	if(ndps_in_scaf > 1) {
		this->cluster_useful_ndps++;
		// Is this scaf part of the cluster after the addition of the dp?
		if(is_scaf_part_of(ndps_in_cluster+1, ndps_in_scaf)) {
			// Is this scaf part of the cluster before the addition of the dp?
			if(!is_scaf_part_of(ndps_in_cluster, ndps_in_scaf)) {
				ndps_in_scafs_that_belong += ndps_in_cluster;
				scafs_in++;
			}
			ndps_in_scafs_that_belong++;
		}
	}

	cluster += dp;

	if(ndps_in_scafs_that_belong < 0 || ndps_in_scafs_that_belong > cluster.ndps()) {
		cerr << "Fatal error (1): illegal number of dps in scafs that belong: " << ndps_in_scafs_that_belong << endl;
		throw exception();
	}

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const ClusterSeparatorSplitScafs::ClusterStats& ClusterSeparatorSplitScafs::ClusterStats::operator -= (size_t dp)
{
	auto	mit = cluster.scafs_find(this->scaf_db.dp2scaf(dp));
	size_t	ndps_in_cluster = (mit != cluster.scafs_end())? mit->second : 0;
	size_t	ndps_in_scaf = scaf_db.ndps(this->scaf_db.dp2scaf(dp));

	if(ndps_in_scaf > 1) {
		this->cluster_useful_ndps--;
		// Is this scaf part of the cluster after the addition of the dp?
		if(is_scaf_part_of(ndps_in_cluster, ndps_in_scaf)) {
			// Is this scaf part of the cluster before the addition of the dp?
			if(!is_scaf_part_of(ndps_in_cluster-1, ndps_in_scaf)) {
				ndps_in_scafs_that_belong -= ndps_in_cluster;
				scafs_in--;
			}
			else
				ndps_in_scafs_that_belong--;
		}
	}
	cluster -= dp;

	if(ndps_in_scafs_that_belong < 0 || ndps_in_scafs_that_belong > cluster.ndps()) {
		cerr << "Fatal error (2): illegal number of dps in scafs that belong: " << ndps_in_scafs_that_belong << endl;
		throw exception();
	}

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ostream& operator << (ostream& os, ClusterSeparatorSplitScafs& cs)
{
        os << " Best separation:\t" << "dimension " << cs.get_separating_dimension() << " ("
                << cs.clustering_db.get_dimension(cs.get_separating_dimension()).get_dimension_name() << ")\t"
                << "Value " << cs.get_separating_value() << "\t"
                << "best split scaf ratio " << cs.get_best_split_scaf_ratio() << endl;

	return os;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterSeparatorSplitScafs::separate_dimension(size_t dimension_index)
{
	semaphore.wait();
	const Dimension& dimension = clustering_db.get_dimension(dimension_index);
	const set<double>& values = dimension.get_values();
	int values_checked = 0;

	// Initialize - curr_cluster1 will have no points while curr_cluster2 will have all dps
	Cluster 					curr_cluster1(scaf_db);
	Cluster 					curr_cluster2(scaf_db);
	ClusterSeparatorSplitScafs::IsScafIn		in_operator(fraction_dps_to_be_considered_in);
	ClusterSeparatorSplitScafs::ClusterStats	curr_cluster_stats1(curr_cluster1, scaf_db, in_operator);
	ClusterSeparatorSplitScafs::ClusterStats	curr_cluster_stats2(curr_cluster2, scaf_db, in_operator);
	vector<Dimension::DataPoint>			dps;

	set<size_t> useful_scafs;
	for(vector<Dimension::DataPoint>::const_iterator vit = dimension.begin(); vit!=dimension.end(); vit++) {
		curr_cluster_stats2 += vit->dp;
		dps.push_back(*vit);
		if(scaf_db.ndps(scaf_db.dp2scaf(vit->dp)) > 1)
			useful_scafs.insert(scaf_db.dp2scaf(vit->dp));
	}

	if(curr_cluster_stats2.num_scafs_in_this_cluster() != useful_scafs.size()) {
			cerr << "Error: # useful scafs=" << useful_scafs.size() << " but in cluster 2: " << curr_cluster_stats2.num_scafs_in_this_cluster() << endl;
		semaphore.notify();
		return;
	}
	sort(dps.begin(), dps.end(), comp_by_value);
	ClusterSeparatorSplitScafs::ClusteringResult	best_dimension_separation;

	for(auto vit=dps.begin(); vit!=dps.end();) {
		double value = vit->value;
		while(vit!=dps.end() && vit->value==value) {
			curr_cluster_stats1 += vit->dp;
			curr_cluster_stats2 -= vit->dp;
			vit++;
		}

		// Now do what we actually need to do in this specific class
		this->separate_specific(curr_cluster_stats1, curr_cluster_stats2, dimension_index, value, useful_scafs.size(), best_dimension_separation);
	}

	// Only a single access to the variables common to all threads
	unique_lock<mutex> lck(update_mtx);
	if(best_dimension_separation.is_better(*best_separation))
		best_separation->copy(best_dimension_separation);

	semaphore.notify();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterSeparatorSplitScafs::separate_specific(ClusterStats& curr_cluster1, ClusterStats& curr_cluster2, size_t dimension_index, double tested_value, size_t nuseful_scafs,
	ClusterSeparatorSplitScafs::ClusteringResult& best_dimension_separation)
{
	// First we will move the dps from cluster2 to cluster1. This includes updating in1, in2, separated, cluster1_useful_ndps and cluster2_useful_ndps.
	size_t separated = nuseful_scafs - curr_cluster1.num_scafs_in_this_cluster() - curr_cluster2.num_scafs_in_this_cluster();

	// Which one is smaller in terms of the number of datapoints that will be considered 
	// for the clustering?
	if((curr_cluster1.ndps_this_cluster_scafs() < this->cluster_ndps_threshold) || (curr_cluster2.ndps_this_cluster_scafs() < this->cluster_ndps_threshold))
		return;

	// Now pick the smaller cluster, we will compute sensitivity and specificity on that one
	const ClusterStats& small_cluster = (curr_cluster1.num_useful_ndps() < curr_cluster2.num_useful_ndps())? curr_cluster1 : curr_cluster2;
	const ClusterStats& large_cluster = (curr_cluster1.num_useful_ndps() >= curr_cluster2.num_useful_ndps())? curr_cluster1 : curr_cluster2;

	double separated_to_cluster1_ratio = double(separated)/double(small_cluster.num_scafs_in_this_cluster());
	double cluster_size_ratio = double(small_cluster.num_scafs_in_this_cluster())/large_cluster.num_scafs_in_this_cluster();
	if(cluster_size_ratio<1)
		cluster_size_ratio = 1/cluster_size_ratio;

	ClusterSeparatorSplitScafs::ClusteringResult	cr(dimension_index, tested_value, separated_to_cluster1_ratio, cluster_size_ratio);
	if(cr.is_better(best_dimension_separation))
		best_dimension_separation.copy(cr);
}
