#ifndef CLUSTERINGDATA_H
#define CLUSTERINGDATA_H

#include <vector>
#include <ostream>
#include <string>
#include <map>
#include <unordered_set>
#include "ScafDpData.h"
#include "Dimension.h"

using namespace std;

typedef vector<double> ValueVector;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterData
// This class stores data that is used for clustering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterData {
public:
	ClusterData(string esom_lrn_file, const ScafDpData& scaf_db);
	ClusterData(const ClusterData& cluster_data, const unordered_set<size_t>& subset_dps);
	~ClusterData();
	size_t			ndimensions() const					{return (dimensions.size()-1);}
	size_t			ndps() const						{return dps.size();}
	const Dimension&	get_dimension(size_t index) const;
	const Dimension&	get_dimension(string dname) const;
	const unordered_set<size_t>&	datapoints() const				{return dps;}
	double			get_value(size_t datapoint, size_t dimension) const	{return dimensions[dimension]->get_value(datapoint);}
protected:
	vector<Dimension*>	dimensions;
	unordered_set<size_t>	dps;
};

#endif 	//CLUSTERINGDATA_H
