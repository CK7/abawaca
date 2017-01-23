#ifndef CLUSTERSEPARATORBYSENSITIVITYSPECIFICITY_H
#define CLUSTERSEPARATORBYSENSITIVITYSPECIFICITY_H

#include "ClusterSeparator.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterSeparatorBySensitivitySpecificity
// Implements the basic strategy that was successfully attempted in the initial PERL version of this algorithm
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparatorBySensitivitySpecificity : public ClusterSeparator {
protected:
	class ClusterStats;
	class ClusteringResult : public ClusterSeparator::ClusteringResult {
	public:
	        ClusteringResult() : ClusterSeparator::ClusteringResult(), specificity(0), sensitivity(0)								{}
	        ClusteringResult(size_t d, double v, double sens, double spec) : ClusterSeparator::ClusteringResult(d, v), specificity(spec), sensitivity(sens)		{}
		void			reset()                                         {sensitivity = specificity = 0; ClusterSeparator::ClusteringResult::reset();}
		void			copy(const ClusterSeparator::ClusteringResult& cr)
	        {
	                this->sensitivity = dynamic_cast<const ClusteringResult*>(&cr)->sensitivity;
	                this->specificity = dynamic_cast<const ClusteringResult*>(&cr)->specificity;
			ClusterSeparator::ClusteringResult::copy(cr);
	        }
	        // Is this one better than cr?
		bool			is_better(const ClusterSeparator::ClusteringResult& cr_) const
			{
				const ClusteringResult& cr = *dynamic_cast<const ClusteringResult*>(&cr_);
				if(sensitivity*specificity != cr.sensitivity*cr.specificity) return (sensitivity*specificity >= cr.sensitivity*cr.specificity);
				return ClusterSeparator::ClusteringResult::is_better(cr);
			}
		bool			is_legal() const
			{
				return ((sensitivity>=sensitivity_threshold) && (specificity>=specificity_threshold) &&
			                (sensitivity+specificity>=sum_threshold) && (sensitivity*specificity>=prodduct_threshold));
			}
		double			get_sensitivity() const				{return sensitivity;}
		double			get_specificity() const				{return specificity;}
		void			print_separation_params(ostream& os) const	{ClusterSeparator::ClusteringResult::print_separation_params(os); os << ", Specificity: " << specificity << ", sensitivity:" << sensitivity;}
	protected:
	        double			sensitivity;
	        double			specificity;
		double			sensitivity_threshold = 0.8;
		double			specificity_threshold = 0.8;
		double			prodduct_threshold = 0.8;
		double			sum_threshold = 1.6;
	};
public:
	ClusterSeparatorBySensitivitySpecificity(const ScafDpData& scaf_db, const SCGdb& scg_db, const ClusterData& clustering_db,  Semaphore& threads_semaphore) :
		ClusterSeparator(scaf_db, scg_db, clustering_db, threads_semaphore, new ClusteringResult())
		{}
	double			get_best_sensitivity() const		{return dynamic_cast<ClusteringResult*>(best_separation)->get_sensitivity();}
	double			get_best_specificity() const		{return dynamic_cast<ClusteringResult*>(best_separation)->get_specificity();}
	friend ostream& 	operator << (ostream& os, ClusterSeparatorBySensitivitySpecificity& cs);
protected:
	void			separate_dimension(size_t dimension_index) override;
	void			separate_specific(ClusterStats& curr_cluster1, ClusterStats& curr_cluster2, size_t dimension_index, double tested_value, ClusteringResult& best_dimension_separation);
//	void 			separate_specific(Cluster& curr_cluster1, Cluster& curr_cluster2, size_t cluster1_useful_ndps, size_t cluster2_useful_ndps, size_t dimension_index, double tested_value);
};

ostream& operator << (ostream& os, ClusterSeparatorBySensitivitySpecificity& cs);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterSeparatorBySensitivitySpecificity::ClusterStats : public ClusterSeparator::ClusterStats {
public:
        ClusterStats(Cluster& _cluster, const ScafDpData& _scaf_db) : ClusterSeparator::ClusterStats(_cluster, _scaf_db) {}
        const ClusterSeparator::ClusterStats&	add_dp(size_t dp) override;
        const ClusterSeparator::ClusterStats&	remove_dp(size_t dp) override;
        const ClusterSeparator::ClusterStats&	assign_scaf(size_t scaf) override;
        const ClusterSeparator::ClusterStats&	unassign_scaf(size_t scaf) override;
	int                     num_useful_ndps() const override		{return (true_positives+false_positives);}
        int                     get_true_positives() const			{return true_positives;}
        int                     get_total_dps_for_assigned_scafs() const	{return total_dps_for_assigned_scafs;}
        int                     get_false_positives() const			{return false_positives;}
	double			get_specificity() const				{return double(true_positives)/(true_positives+false_positives);}
	double			get_sensitivity() const				{return double(true_positives)/total_dps_for_assigned_scafs;}
protected:
        int                     true_positives = 0;	// # of dps in cluster for scafs assigned to cluster
        int                     total_dps_for_assigned_scafs = 0;
	int			false_positives = 0;	// # of dps in cluster for scafs not assigned to this cluster
};

#endif //CLUSTERSEPARATORBYSENSITIVITYSPECIFICITY_H
