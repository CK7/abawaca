/*
 * Sequence.h
 *
 *  Created on: 05/Feb/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include "bio_exceptions.h"

using namespace std;

namespace Bio {

typedef enum {forward=1, reverse=2} Direction;

/***************************************************************************************************************************
 * Bio::Sequence
 * This class is a container for a biological (DNA, RNA, protein, else?) sequence. Real sequence functionality depends on
 * the class T.
 ***************************************************************************************************************************/
template <class T>
class Sequence {
public:
	Sequence(string display_id, string description, T sequence) : m_display_id(display_id), m_description(description),
		m_sequence(sequence) 								{}
	virtual 		~Sequence()							{}
	const string&		display_id() const						{return m_display_id;}
	const string&		desc() const							{return m_description;}
	const T&		seq() const							{return m_sequence;}
	T&			seq()								{return m_sequence;}
	bool			operator < (const Sequence<T>& other) const			{return (this->display_id() < other.display_id());}
protected:
	string	m_display_id;
	string	m_description;
	T	m_sequence;
public:
};

}
#endif /* SEQUENCE_H_ */
