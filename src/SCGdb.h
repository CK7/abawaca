#ifndef SCGDB_H
#define SCGDB_H
#include <string>
#include <map>
#include <set>
#include "ScafDpData.h"

using namespace std;

class SCGdb {
public:
	// _gene2scg_file can be an empty string
	SCGdb(const ScafDpData& _scaf_db, string _gene2scg_file, string _scg_list_file) :gene2scg_file(_gene2scg_file), scg_list_file(_scg_list_file),
		scaf_db(_scaf_db) {this->load_scgs();}
	size_t				num_unique_scgs(const set<size_t>& scafs) const;
	// The following function will put into result of SCG and their counts in the provided set of scaffolds
	void				get_scg2count(const set<size_t>& scafs, map<string, size_t>& result) const;
	double				average_num_copies_for_unique_scgs(const set<size_t>& scafs) const;
	size_t				get_total_num_scgs() const			{return total_num_scgs;}
protected:
	void load_scgs();
protected:
	string				gene2scg_file;
	string				scg_list_file;
	map<size_t, set<string> >	scaf2scgs;
	int				total_num_scgs = 0;
	const ScafDpData&		scaf_db;
};

#endif // SCGDB_H
