#ifndef CLUSTERINGDATA_H
#define CLUSTERINGDATA_H

#include <vector>
#include <ostream>
#include <string>
#include <map>
#include "Dimension.h"

using namespace std;

typedef vector<double> ValueVector;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterData
// This class stores data that is used for clustering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterData {
public:
	ClusterData(const char* esom_lrn_file);
	ClusterData(const ClusterData& cluster_data, const set<size_t>& subset_dps);
	~ClusterData();
	size_t			ndimensions() const					{return (dimensions.size()-1);}
	size_t			ndps() const						{return dps.size();}
	const Dimension&	get_dimension(size_t index) const			{return *(dimensions[index]);}
	const set<size_t>&	datapoints() const					{return dps;}
	double			get_value(size_t datapoint, size_t dimension) const	{return dimensions[dimension]->get_value(datapoint);}
protected:
	vector<Dimension*>	dimensions;
	set<size_t>		dps;
};

#endif 	//CLUSTERINGDATA_H
