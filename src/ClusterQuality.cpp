#include "ClusterQuality.h"
#include <vector>
#include <math.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void standard_deviation(vector<pair<size_t, double> > data, double& mean, double& stdev)
{
	if(data.size() == 0) {
		mean = stdev = -1;
		return;
	}

	size_t	total_size = 0;
	mean=0;
	stdev=0;

	for(auto it=data.begin(); it!=data.end(); it++) {
		mean += it->first*it->second;
		total_size += it->first;
	}
	mean /= total_size;

	for(auto it=data.begin(); it!=data.end(); it++)
		stdev += it->first*(it->second-mean)*(it->second-mean);

	stdev = sqrt(stdev/(total_size-1));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterQuality::dimension(const Cluster& cluster, double& avg, double& stdev, const Dimension& dimension) const
{
	vector<pair<size_t, double> >   d;

	for(auto it=cluster.assigned_scafs_begin(); it!=cluster.assigned_scafs_end(); it++) {
		const set<const ScafDpData::DataPoint*>& scaf_dps = scaf_db.get_scaf(*it)->get_scaf_dps();
		for(auto sit=scaf_dps.begin(); sit!=scaf_dps.end(); sit++) {
			d.push_back(pair<size_t, double>(1, dimension.get_value((*sit)->get_iid())));
		}
	}
	standard_deviation(d, avg, stdev);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterQuality::scg(const Cluster& cluster, double& nunique, double& avg) const
{
	nunique = scg_db.num_unique_scgs(cluster.get_assigned_scafs());
	avg = scg_db.average_num_copies_for_unique_scgs(cluster.get_assigned_scafs());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterQuality::gc(const Cluster& cluster, double& avg_gc, double& gc_stdev) const
{
	vector<pair<size_t, double> >	gc;

	for(auto it=cluster.assigned_scafs_begin(); it!=cluster.assigned_scafs_end(); it++) {
		const ScafDpData::Seq* seq_obj = scaf_db.get_scaf(*it);
		if(seq_obj != NULL)
			gc.push_back(pair<size_t, double>(seq_obj->get_seq_obj()->seq().size(), seq_obj->get_gc()));
	}
	standard_deviation(gc, avg_gc, gc_stdev);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClusterQuality::cvg(const Cluster& cluster, double& avg_cvg, double& cvg_stdev) const
{
	vector<pair<size_t, double> >	cvg;

	for(auto it=cluster.assigned_scafs_begin(); it!=cluster.assigned_scafs_end(); it++) {
		const ScafDpData::Seq* seq_obj = scaf_db.get_scaf(*it);
		if(seq_obj != NULL) {
			cvg.push_back(pair<size_t, double>(seq_obj->get_seq_obj()->seq().size(), seq_obj->get_cvg()));
		}
	}
	standard_deviation(cvg, avg_cvg, cvg_stdev);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t ClusterQuality::total_size(const Cluster& cluster) const
{
	size_t ts = 0;
	for(auto it=cluster.assigned_scafs_begin(); it!=cluster.assigned_scafs_end(); it++) {
		const Bio::DNASequence* seq_obj =  scaf_db.get_scaf(*it)->get_seq_obj();
		if(seq_obj != NULL)
			ts += seq_obj->seq().size();
	}
	return ts;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t ClusterQuality::num_scaffolds(const Cluster& cluster) const
{
	return cluster.get_assigned_scafs().size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClusterQuality::is_split_better(/*const Cluster& cluster, */const Cluster& sub_cluster1, const Cluster& sub_cluster2) const
{
	map<string, size_t>	cluster_scgs, sub_cluster1_scgs, sub_cluster2_scgs;
	size_t	sub_cluster1_size=0, sub_cluster2_size=0;

//	auto cluster_scafs = cluster.get_assigned_scafs();
	auto sub_cluster1_scafs = sub_cluster1.get_assigned_scafs();
	auto sub_cluster2_scafs = sub_cluster2.get_assigned_scafs();
	set<size_t> cluster_scafs(sub_cluster1_scafs);
	cluster_scafs.insert(sub_cluster2_scafs.begin(), sub_cluster2_scafs.end());

	scg_db.get_scg2count(cluster_scafs, cluster_scgs);
	scg_db.get_scg2count(sub_cluster1_scafs, sub_cluster1_scgs);
	scg_db.get_scg2count(sub_cluster2_scafs, sub_cluster2_scgs);

	for(auto it=sub_cluster1_scafs.begin(); it!=sub_cluster1_scafs.end(); it++)
		sub_cluster1_size += scaf_db.get_scaf(*it)->get_seq_size();

	for(auto it=sub_cluster2_scafs.begin(); it!=sub_cluster2_scafs.end(); it++)
		sub_cluster2_size += scaf_db.get_scaf(*it)->get_seq_size();

	// If one of the clusters has no SCGs then decide based on the total cluster size.
	if((sub_cluster1_scgs.size() == 0) || (sub_cluster2_scgs.size() == 0)) {
		return ((sub_cluster1_size >= 500000) && (sub_cluster2_size >= 500000));
	}

	// If the two clusters have complementary sets of SCGs - don't split further.
	set<string> scgs;
	for(auto it=sub_cluster1_scgs.begin(); it!=sub_cluster1_scgs.end(); it++) {
		if(sub_cluster2_scgs.find(it->first) != sub_cluster2_scgs.end())
			scgs.insert(it->first);
	}

	// If the overlap between the two clusters is small then this is probably an over-splitting.
	if(double(scgs.size())/((sub_cluster1_scgs.size()<sub_cluster2_scgs.size())? sub_cluster1_scgs.size() : sub_cluster2_scgs.size()) < 0.2)
		return false;

	// Consider more criteria

	return true;
}
