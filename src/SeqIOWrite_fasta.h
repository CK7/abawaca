/*
 * SeqIO_fasta.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIOWRITE_FASTA_H_
#define SEQIOWRITE_FASTA_H_

#include "SeqIOWrite.h"
#include "common.h"
#include "bio_exceptions.h"
#include <ctype.h>

using namespace std;

namespace Bio {

/****************************************************************************************************************
 * SeqIOWrite_fasta
 * This is an implementation of a fasta file writer. Specific sequence and its type type (sequence, read/DNA, RNA, protein) 
 * is determined by the supplied type.
 ****************************************************************************************************************/
template <class S>
class SeqIOWrite_fasta : public SeqIOWrite<S> {
public:
	SeqIOWrite_fasta(string fasta_file=string("")) : SeqIOWrite<S>(fasta_file), _line_length(60)	{}
	void						open(string file)			{this->SeqIOWrite<S>::open(file); _line_length = 60; }
	void 						write_seq(const S& seq);
	size_t						size_length() const			{return _line_length;}
	void						size_length(size_t l)			{_line_length = l;}
protected:
	size_t	_line_length;
};

/****************************************************************************************************************/
template <class S>
void SeqIOWrite_fasta<S>::write_seq(const S& seq)
{
	this->m_fp << '>' <<  seq.display_id() << ' ' <<  seq.desc() << endl;
	// A conversion from seq.seq() to string is assumed to exist!
	string	str = string(seq.seq());
	this->m_line_num++;
	for(size_t i=0; i<=str.size(); i+=_line_length) {
		this->m_fp << str.substr(i, (i+_line_length<=str.size())? _line_length : str.size()-i+1) << endl;
		this->m_line_num++;
	}
	this->m_seq_num++;
}

}

#endif /* SEQIO_FASTA_H_ */
