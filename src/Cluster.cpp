#include "Cluster.h"

#include <iostream>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cluster::Cluster(const ScafDpData& _scaf_db, const Dimension& dimension, const Op& op, double threshold) : scaf_db(_scaf_db)
{
	for(vector<Dimension::DataPoint>::const_iterator vit=dimension.begin(); vit!=dimension.end(); vit++)
		if(op(vit->value, threshold))
			*this += vit->dp;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t Cluster::ndps(size_t scaf) const 
{
	map<size_t, size_t>::const_iterator mit = scaf2ndps.find(scaf);

	if(mit==scaf2ndps.end())
		return 0;
 
	return mit->second;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const Cluster& Cluster::operator += (size_t dp)
{
	if(dps.find(dp) != dps.end())
		return *this;
	dps.insert(dp);

	size_t scaf = scaf_db.dp2scaf(dp);
	map<size_t, size_t>::iterator mit = scaf2ndps.find(scaf);
	if(mit == scaf2ndps.end())
		scaf2ndps.insert(pair<size_t, size_t>(scaf, 1));
	else
		mit->second++;

	return *this;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const Cluster& Cluster::operator -= (size_t dp)
{
	if(dps.find(dp) == dps.end())
		return *this;
	dps.erase(dp);

	size_t scaf = scaf_db.dp2scaf(dp);
	map<size_t, size_t>::iterator mit = scaf2ndps.find(scaf);
	if(mit->second == 1)
		scaf2ndps.erase(mit);
	else
		mit->second--;

	return *this;
}
