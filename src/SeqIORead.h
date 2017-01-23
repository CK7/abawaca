/*
 * SeqIORead.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIOREAD_H_
#define SEQIOREAD_H_

#include "SeqIO.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;

namespace Bio {

/*****************************************************************************************************************
 * SeqIORead
 * This class imitates BioPerl's Bio::SeqIO with "read" mode, i.e. it only has the next_seq() function but not
 * the write_seq(). Other differences:
 * - The class itself is abstract, user needs to decide which version it picks.
 ****************************************************************************************************************/
template <class S>
class SeqIORead : public SeqIO {
public:
	SeqIORead();
	SeqIORead(string file);
	void						open(string file);
	void						close()					{this->SeqIO::close(); if(m_fp) {fclose(m_fp); m_fp=NULL;}}
	virtual 					~SeqIORead()				{this->close(); free(m_next_line);}
	void						restart();
	// next_seq() returns NULL when no more sequences are left. It is the user's responsibility to destroy
	// the returned object.
	virtual S* 					next_seq() = 0;
protected:
        char*						getline(bool skip_empty=false);
protected:
	FILE*						m_fp;
        char*						m_next_line;
        size_t						m_next_line_size;
};

/****************************************************************************************************************/
template<class S>
SeqIORead<S>::SeqIORead() : SeqIO(), m_fp(NULL), m_next_line((char*)calloc(1024, sizeof(char))), 
	m_next_line_size(1024)
{
}

/****************************************************************************************************************/
template<class S>
SeqIORead<S>::SeqIORead(string file) : SeqIO(file), m_fp(NULL), m_next_line((char*)calloc(1024, sizeof(char))), 
	m_next_line_size(1024)
{
	this->open(file);
}

/****************************************************************************************************************/
template<class S>
void SeqIORead<S>::open(string file)
{
	this->SeqIO::open(file);
	if(!(m_fp = fopen(this->file_name().c_str(), "r")))
		throw Bad_file(file, "Failed to open file", __PRETTY_FUNCTION__);
} 

/****************************************************************************************************************/
template<class S>
void SeqIORead<S>::restart()
{
	SeqIO::restart();

	if(m_fp)
		fclose(m_fp);

	if(!(m_fp = fopen(this->file_name().c_str(), "r")))
		throw Bio::Bad_file(this->file_name(), "Failed to re-open file", __PRETTY_FUNCTION__);
}

/****************************************************************************************************************/
template <class S>
char* SeqIORead<S>::getline(bool skip_empty/*=false*/)
{
	if(feof(this->m_fp))
		return NULL;

	do {
		this->m_next_line[0] = 0;
		// We do the following because we don't know whether the size of next_line is sufficient for storing the next line.
		// If not we will increase the size of next_line until we have enough space. Eventually next_line will be able to
		// store the longest line in the file.
		char    *p = this->m_next_line;
		size_t  s = this->m_next_line_size;
		while(fgets(p, s, this->m_fp)) {
			if((*(p+strlen(p)-1) == '\n') || feof(this->m_fp))
				break;
			if(strlen(p) != s-1) {
				// This is weird!
			}
			// Else - we did not read the whole line, need to increase next_line's size and continue reading
			this->m_next_line_size *= 2;
			this->m_next_line = (char*)realloc(this->m_next_line, this->m_next_line_size);
			p = this->m_next_line+strlen(this->m_next_line);
			s = this->m_next_line_size - strlen(this->m_next_line);
		}
		p = this->m_next_line;
		while(*p && isspace(*p))
			p++;
		// If we reached the end of the line and everything we saw is just spaces...
		if(!(*p) && skip_empty)
			this->m_next_line[0] = 0;
		this->m_line_num++;
	} while((this->m_next_line[0] == 0) && !feof(this->m_fp));

	if(m_next_line[0] == 0)
		return NULL;
	return m_next_line;
}

}

#endif /* SEQIOREAD_H_ */
