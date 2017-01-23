/*
 * SeqIOWrite.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIOWRITE_H_
#define SEQIOWRITE_H_

#include "SeqIO.h"

using namespace std;

namespace Bio {

/*****************************************************************************************************************
 * SeqIOWrite
 * This class imitates BioPerl's Bio::SeqIO with "write" mode, i.e. it only has the write_seq() function but not
 * the next_seq(). Other differences:
 * - The class itself is abstract, user needs to decide which version it picks.
 ****************************************************************************************************************/
template <class S>
class SeqIOWrite : public SeqIO {
public:
	SeqIOWrite(string file);
	virtual 					~SeqIOWrite()				{this->close();}
	void						close()					{this->SeqIO::close(); if(m_fp.is_open()) m_fp.close();}
	void						open(string file);
	virtual void 					write_seq(const S& seq) = 0;
	void						restart();
protected:
	ofstream					m_fp;
};

/****************************************************************************************************************/
template<class S>
void SeqIOWrite<S>::restart()
{
	SeqIO::restart();

	m_fp.close();
	m_fp.open(this->file_name().c_str());

	if(m_fp.fail())
		throw Bio::Bad_file(this->file_name(), "Failed to re-open file", __PRETTY_FUNCTION__);
}


/****************************************************************************************************************/
template<class S>
SeqIOWrite<S>::SeqIOWrite(string file) : SeqIO(file)
{
	if(file.size() > 0) {
		this->open(file);
	}
}

/****************************************************************************************************************/
template<class S>
void SeqIOWrite<S>::open(string file)
{
	this->close();
	this->SeqIO::open(file);
	m_fp.open(file.c_str());

	if(m_fp.fail())
		throw Bad_file(file, "Failed to open file", __PRETTY_FUNCTION__);
}

}

#endif /* SEQIOWRITE_H_ */
