#ifndef SCAFDPDATA_H
#define	SCAFDPDATA_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "String.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ScafDpData {
public:
	class DataPoint;
	class Seq;
public:
	ScafDpData(string esom_names_file, string fasta_file, string info_file);
	~ScafDpData()							{}
	size_t			nscafs() const				{return scafs.size()-1;}
	size_t			ndps() const				{return dps.size()-1;}
	size_t			ndps(size_t scaf_id) const;
	size_t			dp2scaf(size_t dp) const;
	string			scaf_id2name(size_t scaf) const;
	size_t			scaf_name2id(string scaf) const		{if(scaf_name2id_map.find(scaf)!=scaf_name2id_map.end()) return scaf_name2id_map.find(scaf)->second; return 0;}
	size_t			dp_name2id(size_t dp) const		{auto it=dp_name2id_map.find(dp); return (it != dp_name2id_map.end())? it->second : 0;}
	const Seq*		get_scaf(size_t scaf) const		{return (scaf>=scafs.size())? NULL : &(scafs[scaf]);}
	const Seq*		get_scaf(const string& scaf) const	{auto it=scaf_name2id_map.find(scaf); return (it==scaf_name2id_map.end())? NULL : this->get_scaf(it->second);}
//	const DataPoint*	get_dp(size_t id) const			{return (id < dps.size())? &(dps[id]) : NULL;}
	const DataPoint*	get_dp_by_name(size_t dp) const		{auto it=dp_name2id_map.find(dp); return (it != dp_name2id_map.end())? &(dps[it->second]) : NULL;}
protected:
	map<string, size_t>	scaf_name2id_map;
	map<size_t, size_t>	dp_name2id_map;
	vector<DataPoint>	dps;
	vector<Seq>		scafs;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class represents a sequence in the database and contains all its related information
class ScafDpData::Seq {
public:
	Seq() : seq_id(0), seq_obj(NULL), cvg(0)										{}
	Seq(size_t _seq_id, const Bio::DNASequence* _seq_obj, double _cvg=0) : seq_id(_seq_id), seq_obj(_seq_obj), cvg(_cvg)	{}
	const Seq&			operator = (const Seq& s) 								{seq_id=s.seq_id; cvg=s.cvg; seq_obj=s.seq_obj; scaf_dps.erase(scaf_dps.begin(), scaf_dps.end()); scaf_dps.insert(s.scaf_dps.begin(), s.scaf_dps.end()); return *this;}
	const set<const DataPoint*>&	get_scaf_dps() const									{return scaf_dps;}
	const Bio::DNASequence*		get_seq_obj() const									{return seq_obj;}
	double				get_cvg() const										{return cvg;}
	double				get_gc() const										{return Bio::gc(seq_obj->seq());}
	size_t				get_Ns() const										{return Bio::Ns(seq_obj->seq());}
	size_t				get_seq_size() const									{return seq_obj->seq().size();}
protected:
	// dps are added after the creation of the object by ScafDpData
	friend class ScafDpData;
	const Seq&			operator +=(const DataPoint& dp)							{scaf_dps.insert(&dp); return *this;}
protected:
	size_t				seq_id;
	set<const DataPoint*>		scaf_dps;
	const Bio::DNASequence* 	seq_obj;
	double				cvg;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ScafDpData::DataPoint {
public:
	DataPoint(size_t _dp_id, size_t _dp_name, size_t _scaf_id, size_t _start, size_t _end, size_t _Ns) : dp_id(_dp_id), dp_name(_dp_name), scaf_id(_scaf_id), start(_start), end(_end), Ns(_Ns) {}
	DataPoint() : dp_id(0), dp_name(0), scaf_id(0), start(-1), end(-1), Ns(-1) 	{}
	DataPoint&		operator = (const DataPoint& o)=default;
	// That's the internal ID used throughout the program
	size_t			get_iid() const				{return dp_id;}
	// That's the name given to the dp in the .lrn file
	size_t			get_name() const			{return dp_name;}
	size_t			get_scaf_id() const			{return scaf_id;}
	size_t			get_start() const			{return start;}
	size_t			get_end() const				{return end;}
	size_t			get_Ns() const				{return Ns;}
	size_t			length() const				{return (end-start+1);}
	size_t			non_Ns() const				{return (end-start+1-Ns);}
	bool			operator < (const DataPoint& dp) const	{return (this->dp_name < dp.dp_name);}
protected:
	size_t			dp_id;
	size_t			dp_name;
	size_t			scaf_id;
	size_t			start;
	size_t			end;
	size_t			Ns;
};

#endif	// SCAFDPDATA_H
