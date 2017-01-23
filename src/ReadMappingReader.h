/*
 * ReadMappingReader.h
 *
 *  Created on: 05/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 * This version is implemented using low-level C features (FILE*, pointers instead of ifstream and string) in order to 
 * increase speed. As a result this version takes ~40% less time than the first one.
 */

#ifndef READMAPPINGREADER_H_
#define READMAPPINGREADER_H_

#include "common.h"
#include "ReadMapping.h"
#include "SeqIO.h"
#include "bio_exceptions.h"
#include <stdio.h>
#include <vector>
#include <stdlib.h>

using namespace std;

namespace Bio {

/***************************************************************************************************************************
 * ReadMappingReader
 * Base class for classes reading read mapping.
 ***************************************************************************************************************************/
class ReadMappingReader {
public:
	ReadMappingReader()		{}
	virtual 			~ReadMappingReader()		{}
	virtual	ReadMapping*		next_mapping() = 0;
	virtual void			restart() = 0;
	virtual bool			good() const = 0;
	virtual double			percent_entries_read() const = 0;
};

/***************************************************************************************************************************
 * FileReader
 * Base class for classes reading read mapping files of all formats
 ***************************************************************************************************************************/
class FileReader : public ReadMappingReader {
public:
	FileReader(string file);
	virtual 			~FileReader()			{if(m_fp!=NULL) {fclose(m_fp);}}
	void				restart();
	const string&			file_name() const		{return m_file_name;}
	bool				good() const			{return !feof(m_fp);} 
	double				percent_entries_read() const	{if(m_total_to_read==0) return 100.0; return (100000*m_total_read/m_total_to_read)/1000.0;}
protected:
	FILE*				m_fp;
	string				m_file_name;
	uint64_t			m_total_to_read;
	uint64_t			m_total_read;
};

/***************************************************************************************************************************
 * SAMReader
 * Implementation for the SAM format
 ***************************************************************************************************************************/
class SAMReader : public FileReader {
public:
	SAMReader(string sam_file) : FileReader(sam_file)		{this->read_header();}
	ReadMapping*			next_mapping();
	// Currently does nothing, have to decide what to do
	void				read_header();
};

/***************************************************************************************************************************
 * VectorReader
 * Takes ReadMapping from a vector provided by the user. Useful for cases in which a program generates the mapping during
 * runtime (and not through external mapper).
 ***************************************************************************************************************************/
class VectorReader : public ReadMappingReader {
public:
	VectorReader(const vector<const ReadMapping*>& _read_mappings) : ReadMappingReader(), read_mappings(_read_mappings), it(read_mappings.begin())		{}
	ReadMapping*				next_mapping(); //		{return (it==read_mappings.end())? NULL : new ReadMapping(*(*(it++)));}
	void					restart()			{it = read_mappings.begin();}
	bool					good() const			{return (it != read_mappings.end());}
	double					percent_entries_read() const	{if(it==read_mappings.end()) return 100.0; return (100000*(it-read_mappings.begin())/read_mappings.size())/1000.0;}
protected:
	vector<const ReadMapping*> 		read_mappings;
	vector<const ReadMapping*>::const_iterator	it;
};

}
#endif /* READMAPPINGREADER_H_ */
