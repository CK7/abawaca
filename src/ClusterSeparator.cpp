#include "ClusterSeparator.h"
#include <iostream>
#include <algorithm>
#include "misc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClusterSeparatorBySensitivitySpecificity::legal(double sensitivity, double specificity) const
{
	return ((sensitivity>=sensitivity_threshold) && (specificity>=specificity_threshold) && 
		(sensitivity+specificity>=sum_threshold) && (sensitivity*specificity>=prodduct_threshold));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClusterSeparatorBySensitivitySpecificity::best(double sensitivity, double specificity) const
{
	return (sensitivity*specificity>=best_sensitivity*best_specificity);
}

bool comp_by_value (const Dimension::DataPoint& i, const Dimension::DataPoint& j) { return (i.value<j.value); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClusterSeparator::separate()
{
	// Start by cleaning everything from previous runs
	if(cluster1 != NULL)
		delete(cluster1);
	if(cluster2 != NULL)
		delete(cluster2);
	cluster1 = cluster2 = NULL;
	separating_dimension = -1;
	separating_value = -1;
	raw_dps_cluster1.clear();
	raw_dps_cluster2.clear();

	// Do algorithm-specific adjustments
	this->initialize_specific();

	cerr << '[' << get_time() << ']' << " Going over " << clustering_db.ndimensions() << " dimensions" << endl;

	// Go over all dimensions (outern loop)....
	for(size_t dimension_index=1; dimension_index<=clustering_db.ndimensions(); dimension_index++) {
		const Dimension& dimension = clustering_db.get_dimension(dimension_index);
		const set<double>& values = dimension.get_values();
		cerr << '[' << get_time() << ']' << " Dimension " << dimension_index << "/" << clustering_db.ndimensions() << " (" << clustering_db.get_dimension(dimension_index).get_dimension_name() << "), " << values.size() << " values to check" << endl;
		int values_checked = 0;

		// Initialize - curr_cluster1 will have no points while curr_cluster2 will have all dps
		Cluster curr_cluster1(scaf_db);
		Cluster curr_cluster2(scaf_db);
		vector<Dimension::DataPoint>	dps;
		size_t cluster1_useful_ndps = 0, cluster2_useful_ndps = 0;

		for(vector<Dimension::DataPoint>::const_iterator vit = dimension.begin(); vit!=dimension.end(); vit++) {
			if(scaf_db.ndps(scaf_db.dp2scaf(vit->dp)) > 1)
				cluster2_useful_ndps++;
			curr_cluster2 += vit->dp;
			dps.push_back(*vit);
		}
		sort(dps.begin(), dps.end(), comp_by_value);

		for(vector<Dimension::DataPoint>::const_iterator vit=dps.begin(); vit!=dps.end();) {
			double value = vit->value;
			while(vit!=dps.end() && vit->value==value) {
				if((vit-dps.begin()+1)%2000 == 0)
					cerr << '[' << get_time() << ']' << " " << (vit-dps.begin()+1)  << endl;
				curr_cluster1 += vit->dp;
				curr_cluster2 -= vit->dp;
				if(scaf_db.ndps(scaf_db.dp2scaf(vit->dp)) > 1) {
					cluster1_useful_ndps++;
					cluster2_useful_ndps--;
				}
				vit++;
			}

			// Now do what we actually need to do in this specific class
			this->separate_specific(curr_cluster1, curr_cluster2, cluster1_useful_ndps, cluster2_useful_ndps, dimension_index, value);
		}
	}
	// If no significant separation was found simply return false
	if(cluster1 == NULL)
		return false;

	// Before switching datapoints between the clusters keep the original sets
	for(set<size_t>::const_iterator sit = cluster1->dps_begin(); sit != cluster1->dps_end(); sit++)
		raw_dps_cluster1.insert(*sit);
	for(set<size_t>::const_iterator sit = cluster2->dps_begin(); sit != cluster2->dps_end(); sit++)
		raw_dps_cluster2.insert(*sit);

	// Now decide, for each scaffold, to which cluster it belong. Cluster1 will contain all scaffolds who have at least half
	// of their datapoints in it, cluster2 will contain the rest
	for(map<size_t, size_t>::const_iterator mit = cluster1->scafs_begin(); mit != cluster1->scafs_end(); mit++)
		if(2*mit->second >= scaf_db.ndps(mit->first))
			cluster1->add_assigned_scaf(mit->first);
	for(map<size_t, size_t>::const_iterator mit = cluster2->scafs_begin(); mit != cluster2->scafs_end(); mit++)
		if(2*mit->second > scaf_db.ndps(mit->first))
			cluster2->add_assigned_scaf(mit->first);

	// Now move the datapoints based on their scaffold assignment
	set<size_t> dps;
	for(set<size_t>::const_iterator sit = cluster1->dps_begin(); sit != cluster1->dps_end(); sit++)
		if(!cluster1->assigned(scaf_db.dp2scaf(*sit)))
			dps.insert(*sit);

	for(set<size_t>::const_iterator sit = dps.begin(); sit != dps.end(); sit++) {
		*cluster1 -= *sit;
		*cluster2 += *sit;
	}

	dps.clear();
	for(set<size_t>::const_iterator sit = cluster2->dps_begin(); sit != cluster2->dps_end(); sit++)
		if(!cluster2->assigned(scaf_db.dp2scaf(*sit)))
			dps.insert(*sit);

	for(set<size_t>::const_iterator sit = dps.begin(); sit != dps.end(); sit++) {
		*cluster1 += *sit;
		*cluster2 -= *sit;
	}

	// Do algorithm-specific post processing and return the reply 
	return this->post_specific();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterSeparatorBySensitivitySpecificity::separate_specific(Cluster& curr_cluster1, Cluster& curr_cluster2, size_t cluster1_useful_ndps, size_t cluster2_useful_ndps,
			size_t dimension_index, double tested_value)
{
	// Which one is smaller in terms of the number of datapoints that will be considered 
	// for the clustering?
	if((cluster1_useful_ndps < 40) || (cluster2_useful_ndps < 40))
		return;
			
	// Now pick the smaller cluster, we will compute sensitivity and specificity on that one
	const Cluster& small_cluster = (cluster1_useful_ndps < cluster2_useful_ndps)? curr_cluster1 : curr_cluster2;

	// Go over all dps from the small cluster, decide for each scaffold whether it belongs to this cluster or not,
	// and count true/false positives/negatives 
	size_t true_positives=0, false_positives=0, true_negatives=0;
	map<size_t, size_t>::const_iterator mit, end;
	for(mit=small_cluster.scafs_begin(), end=small_cluster.scafs_end(); mit!=end; mit++) {
		size_t scaf_total_ndps = scaf_db.ndps(mit->first);
		if(scaf_total_ndps < 2)
			continue;
		if(2*mit->second >= scaf_total_ndps) {
			true_positives += mit->second;
			true_negatives += (scaf_total_ndps-mit->second);
		}
		else {
			false_positives += mit->second;
		}
	}

	// Compute specificity and sensitivity and check if this one is the best separation so far
	double sensitivity=double(true_positives)/(true_positives+true_negatives);
	double specificity=double(true_positives)/(true_positives+false_positives);
	if(this->best(sensitivity, specificity)) {
		// If this combination of specificity and sensitivity is the best we would like to keep its
		// parameters even if the clusters are illegal
		best_sensitivity = sensitivity; 
		best_specificity = specificity;
		separating_dimension = dimension_index;
		separating_value = tested_value;

		// ... but we would like to keep the clusters only if the splitting is legal
		if(this->legal(sensitivity, specificity)) {
			if(cluster1 != NULL)
				delete(cluster1);
			if(cluster2 != NULL)
				delete(cluster2);
			cluster1 = new Cluster(curr_cluster1);
			cluster2 = new Cluster(curr_cluster2);
			// We want cluster1 to be the smaller and cluster2 to be the larger
			if(cluster1_useful_ndps > cluster2_useful_ndps) {
				Cluster* temp = cluster1;
				cluster1 = cluster2;
				cluster2 = temp;
			}
		}
	}
}
