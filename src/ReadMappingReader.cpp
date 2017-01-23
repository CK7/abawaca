/*
 * ReadMappingReader.cpp
 *
 *  Created on: 05/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include "ReadMappingReader.h"
#include "SeqIO.h"
#include "bio_exceptions.h"
#include <stdlib.h>

namespace Bio {

/***************************************************************************************************************************/
ReadMapping* VectorReader::next_mapping()
{
        ReadMapping* r = NULL;

        if(it==read_mappings.end())
                return NULL;

        try {
                r = new ReadMapping(*(*(it)));
        }
        catch(exception& e) {
		throw Software_bug("Failed to allocate new ReadMapping", "", __PRETTY_FUNCTION__);
        }
	it++;

	return r;
}

/***************************************************************************************************************************/
FileReader::FileReader(string file) : m_fp(NULL), m_file_name(file), m_total_to_read(0),
	m_total_read(0)
{
	if((m_fp = fopen(file.c_str(), "r")) == NULL)
		throw Bio::Bad_file(file, "Failed to open SAM file", __PRETTY_FUNCTION__);
	m_total_to_read = filesize(m_file_name.c_str());
}

/****************************************************************************************************************/
void FileReader::restart()
{
	if(m_fp)
		fclose(m_fp);
	if((m_fp = fopen(this->file_name().c_str(), "r")) == NULL)
		throw Bio::Bad_file(this->file_name(), "Failed to re-open file", __PRETTY_FUNCTION__);
	m_total_to_read = filesize(m_file_name.c_str());
	m_total_read = 0;
}

/***************************************************************************************************************************
 * SAMReader
 ***************************************************************************************************************************/
void SAMReader::read_header()
{
	// At the moment I don't see any reason to keep any of the information provided in the header, i.e. the name of the
	// reference sequences and their lengths
	// Examples:
	// First line -
	// @HD	VN:1.0	SO:unsorted
	// Next lines - @SQ	SN:<ref-name>	LN:<its length>
	// @SQ	SN:310666	LN:2943
/*	m_total_read = 0;
	m_total_to_read = filesize(file.c_str());

	string	next_line;
	ifstream fp;
	fp.open(this->file_name().c_str());

	do {
		getline(m_fp, next_line);
		m_total_to_read -= next_line.size();
	} while(m_fp.good() && ((next_line.size() == 0) || (next_line.at(0) == '@')));*/
}

/***************************************************************************************************************************/
ReadMapping* SAMReader::next_mapping()
{
	// The following loop should have an effect just once - at the beginning. Otherwise it should result with an entry 
	// every time it reads a line 
	static char* 	next_line = (char*)calloc(1024, sizeof(char));
	static size_t	next_line_size = 1024;
	next_line[0] = 0;
	do {
		// We do the following because we don't know whether the size of next_line is sufficient for storing the next line. 
		// If not we will increase the size of next_line until we have enough space. Eventually next_line will be able to 
		// store the longest line in the file.
		char 	*p = next_line;
		size_t	s = next_line_size;
		while(fgets(p, s, m_fp)) {
			if((*(p+strlen(p)-1) == '\n') || !this->good())
				break;
			if(strlen(p) != s-1) {
				// This is weird!
			}
			// Else - we did not read the whole line, need to increase next_line's size and continue reading
			next_line_size *= 2;
			next_line = (char*)realloc(next_line, next_line_size);
			p = next_line+strlen(next_line); 
			s = next_line_size - strlen(next_line);
		}
		if(next_line[0] == '@')
			m_total_to_read -= strlen(next_line);		
	} while(this->good() && ((next_line[0] == 0) || (next_line[0] == '@')));
	if(next_line[0] == 0) {
		return NULL;
	}
//		return new ReadMapping(SAMReader::_end);
	m_total_read += strlen(next_line);
	chomp(next_line);

	return new ReadMapping(next_line);
}

}
