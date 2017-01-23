/*
 * SeqIO.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIO_H_
#define SEQIO_H_

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <ctype.h>
#include "common.h"
#include "bio_exceptions.h"
#include "Sequence.h"

using namespace std;

namespace Bio {

/*****************************************************************************************************************
 * SeqIO
 * This is the base class for SeqIORead and SeqIOWrite.
 ****************************************************************************************************************/
class SeqIO {
public:
	SeqIO(string file=string("")) : m_file_name(file), m_line_num(1), m_seq_num(0)		{}
	virtual ~SeqIO()								{}
	virtual void					close()				{m_file_name=string(""); m_line_num=1; m_seq_num=0;}
	virtual void					open(string file)		{this->close(); m_file_name=file;}	
	string						file_name() const		{return m_file_name;}
	size_t						line_num() const		{return m_line_num;}
	size_t						seq_num() const			{return m_seq_num;}
	virtual void					restart()			{m_seq_num = 0; m_line_num = 1;};
protected:
	string						m_file_name;
	size_t						m_line_num;
	size_t						m_seq_num;
};

}

#endif /* SEQIO_H_ */
