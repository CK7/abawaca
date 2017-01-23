#include <stdlib.h>
#include <iostream>
#include "Dimension.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Dimension::Dimension(const vector<size_t>& dp_ids, const vector<double>& dp_values, size_t _dimension_index, 
	string _dimension_name) : dimension_index(_dimension_index), dimension_name(_dimension_name)
{
	if(dp_values.size() != dp_ids.size()) {
		cerr <<  __FILE__ << " (" << __LINE__ << "): fatal error, number of datapoints (" << dp_ids.size() << ") and values (" 
			<< dp_values.size() << ")for dimension " << dimension_name << " is different" << endl << endl;
		exit(-1);
	}

	vector<size_t>::const_iterator vits = dp_ids.begin();
	vector<double>::const_iterator vitd = dp_values.begin();
	while(vits != dp_ids.end())
		this->insert(*(vits++), *(vitd++));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Dimension::Dimension(const Dimension& dimension, const unordered_set<size_t>& dp_ids) : dimension_index(dimension.dimension_index), 
	dimension_name(dimension.dimension_name)
{
	for(auto vits = dp_ids.begin(); vits != dp_ids.end(); vits++)
		this->insert(*vits, dimension.get_value(*vits));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Dimension::Dimension(size_t _dimension_index, string _dimension_name) : dimension_index(_dimension_index), 
	dimension_name(_dimension_name)
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const Dimension& Dimension::insert(size_t dp, double value)
{
	values.insert(value);
	dp_index.insert(pair<size_t, size_t>(dp, dps.size()));
	dps.push_back(DataPoint(dp, value));

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Dimension::get_value(size_t dp) const
{
	map<size_t, size_t>::const_iterator mit = dp_index.find(dp);
	if(mit == dp_index.end())
		return -1;

	return dps[mit->second].value;
}
