#ifndef DIMENSION_H
#define DIMENSION_H

#include <vector>
#include <set>
#include <string>
#include <map>

using namespace std;

class ClusterData;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dimension
// This class represents a single dimension in a dataset
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Dimension {
public:
	struct DataPoint{
		size_t		dp;
		double		value;
		DataPoint() : dp(-1), value(-1)	{}
		DataPoint(size_t _dp, double _value) : dp(_dp), value(_value)			{}
		bool 				operator < (const DataPoint& other) const	{return (this->dp < other.dp);}
		bool 				operator > (const DataPoint& other) const	{return (this->dp > other.dp);}
		bool 				operator == (const DataPoint& other) const	{return (this->dp == other.dp);}
	};
public:
	size_t					get_dimension_index() const			{return dimension_index;}
	string					get_dimension_name() const			{return dimension_name;}
	// Returns the vector of all values that are represented in the dimension
	const set<double>&			get_values() const				{return values;}
	// Will return the value for dp, or -1 if dp is not found
	double					get_value(size_t dp) const;
	vector<DataPoint>::const_iterator	begin() const					{return dps.begin();}
	vector<DataPoint>::const_iterator	end() const					{return dps.end();}
protected:
	Dimension(const vector<size_t>& dp_ids, const vector<double>& dp_values, size_t dimension_index, string dimension_name);
	Dimension(const Dimension& dimension, const set<size_t>& dp_ids);
	Dimension(size_t dimension_index, string dimension_name);
	const Dimension& insert(size_t dp, double value);
	friend class ClusterData;
protected:
	set<double>			values;
	map<size_t, size_t>		dp_index;
	vector<DataPoint>		dps;
	size_t				dimension_index;
	string				dimension_name;
};

#endif // DIMENSION_H
