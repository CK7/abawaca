/*
 * SeqIORead_fasta.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIOREAD_FASTA_H_
#define SEQIOREAD_FASTA_H_

#include "SeqIORead.h"
#include <ctype.h>

using namespace std;

namespace Bio {

/****************************************************************************************************************
 * SeqIORead_fasta
 * This is an implementation of a fasta file reader. Specific sequence and its type (sequence, read/DNA, RNA, protein) 
 * is determined by the supplied type.
 ****************************************************************************************************************/
template <class S>
class SeqIORead_fasta : public SeqIORead<S> {
public:
	SeqIORead_fasta(string fasta_file);
	S*		 				next_seq();
	void						restart();
};

/****************************************************************************************************************/
template <class S>
SeqIORead_fasta<S>::SeqIORead_fasta(string fasta_file) :
	SeqIORead<S>(fasta_file)
{
	this->getline(true);
}

/****************************************************************************************************************/
template <class S>
void SeqIORead_fasta<S>::restart()
{
	SeqIORead<S>::restart();

	this->getline(true);
}

/****************************************************************************************************************/
template <class S>
S* SeqIORead_fasta<S>::next_seq()
{
	// If reached EOF - return
	if(this->m_next_line[0] == 0)
		return NULL;

	char	*pdisplay_id = this->m_next_line;
	while(*pdisplay_id && isspace(*pdisplay_id))
		pdisplay_id++;
	if((*pdisplay_id != '>') || (*(pdisplay_id+1) == 0) || isspace(*(pdisplay_id+1)))
	{
		char msg[512];
		sprintf(msg, "Was expeecting a header line for the next sequence in a fasta file but got something else, line %u", (unsigned int)this->m_line_num);
		throw Bio::Bad_file(this->file_name(), msg, __PRETTY_FUNCTION__);
	}
	pdisplay_id++;
	char *pdesc = pdisplay_id+1;

	while(*pdesc && !isspace(*pdesc))
		pdesc++;
	if(*pdesc) {
		// Set the end of the sequence name
		*(pdesc++) = 0;
		// Skip the spaces to the beginning of the description (if exists)
		while(*pdesc && isspace(*pdesc))
			pdesc++;
	}

	string	display_id(pdisplay_id), description(pdesc), sequence;
	
	// Now we are ready to read the sequence itself 
	while(this->getline(true)) {
		// Trim before and after
		char *p = this->m_next_line+strlen(this->m_next_line);	// Get p to the \0 at the end of the string
		while(isspace(*(p-1)))
			p--;
		*p = 0;
		p = this->m_next_line;
		while(isspace(*p))
			p++;
		// If this is true - we reached the next entry
		if(*p == '>')
			break;
		sequence += string(p);
	}

	if(sequence.size() == 0) {
		cerr << "Warning: sequence " << display_id << " is empty" << endl;
	}

	this->m_seq_num++;
	return new S(display_id, description, sequence);
}

}

#endif /* SEQIO_FASTA_H_ */
