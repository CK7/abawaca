#ifndef CLUSTER_H
#define CLUSTER_H

#include <map>
#include <set>
#include <unordered_set>
#include "ScafDpData.h"
#include "Dimension.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cluster
// This class represents a single cluster
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Cluster {
public:
	struct Op {virtual bool operator () (double value, double threshold) const = 0;};
	struct Lt : public Op {bool operator () (double value, double threshold) const {return (value < threshold);}};
	struct Le : public Op {bool operator () (double value, double threshold) const {return (value <= threshold);}};
	struct Gt : public Op {bool operator () (double value, double threshold) const {return (value > threshold);}};
	struct Ge : public Op {bool operator () (double value, double threshold) const {return (value >= threshold);}};
public:
	Cluster(const ScafDpData& _scaf_db, const Dimension& dimension, const Op& op, double threshold);
	Cluster(const ScafDpData& _scaf_db) : scaf_db(_scaf_db)				{}

	/* Manipulate datapoints */
	const Cluster&	operator += (size_t dp);
	const Cluster&	operator -= (size_t dp);
	// Iterate through datapoints that belong to this cluster
	unordered_set<size_t>::const_iterator	dps_begin() const			{return dps.begin();}
	unordered_set<size_t>::const_iterator	dps_end() const				{return dps.end();}
	const unordered_set<size_t>&		get_dps() const				{return dps;}
	size_t					ndps() const				{return dps.size();}

	/* Manipulate assigned scafs, i.e. scaffolds that were assigned to the cluster by the user */
	const Cluster&				add_assigned_scaf(size_t scaf)		{assigned_scafs.insert(scaf); return *this;}
	const Cluster&				remove_assigned_scaf(size_t scaf)	{assigned_scafs.erase(scaf); return *this;}
	// Iterate through scaffolds that were assigned to this cluster by the user
	// The following three functions are redundant given that the user can get the whole set... Remove in future versions.
	set<size_t>::const_iterator		assigned_scafs_begin() const		{return assigned_scafs.begin();}
	set<size_t>::const_iterator		assigned_scafs_end() const		{return assigned_scafs.end();}
	bool					assigned(size_t scaf) const		{return (assigned_scafs.find(scaf)!=assigned_scafs.end());}
	const set<size_t>&			get_assigned_scafs() const		{return assigned_scafs;}
	size_t					nassigned_scafs() const			{return assigned_scafs.size();}
	/* Manipulate represented scaffolds, i.e. scaffolds that are represnted by one or more datapoint(s) */
	// Iterate through scaffolds that are represented in this cluster by datapoints
	map<size_t, size_t>::const_iterator	scafs_begin() const			{return scaf2ndps.begin();}
	map<size_t, size_t>::const_iterator	scafs_end() const			{return scaf2ndps.end();}
	map<size_t, size_t>::const_iterator	scafs_find(size_t scaf) const		{return scaf2ndps.find(scaf);}
	size_t					nscafs() const				{return scaf2ndps.size();}
	// Return the number of datapoints that represent scaf in this cluster
	size_t					ndps(size_t scaf) const;
protected:
	const ScafDpData&			scaf_db;
	map<size_t, size_t>			scaf2ndps;
	unordered_set<size_t>			dps;
	set<size_t>				assigned_scafs;
};

#endif // CLUSTER_H
