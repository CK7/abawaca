/*
 * String.cpp
 *
 *  Created on: 05/Feb/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#include "String.h"

namespace Bio {

/**************************************************************************************************************************/
DNAString::DNAString(const string& dna_string)
{
	m_str = dna_string;
/*	for(string::iterator it=m_str.begin(); it!=m_str.end(); it++)
	{
		*it = toupper(*it);
	}*/
	this->validate_string();
}

/**************************************************************************************************************************/
const DNAString& DNAString::operator = (const string& s)
{
	m_str = s;
/*	for(string::iterator it=m_str.begin(); it!=m_str.end(); it++)
	{
		*it = toupper(*it);
	}*/
	this->validate_string();

	return *this;
}

/**************************************************************************************************************************/
void DNAString::validate_string()
{
	const char *p, *pstart = m_str.c_str();
	for(p=pstart; *p; p++) {
		if((*p == 'A') || (*p == 'C') || (*p == 'G') || (*p == 'T') || (*p == 'N')) {
			continue;
		}

		m_str[p-pstart] = toupper(m_str[p-pstart]);

		if((*p != 'A') && (*p != 'C') && (*p != 'G') && (*p != 'T') && (*p == 'N')) {
			throw Illegal_DNAString(m_str, __PRETTY_FUNCTION__, "");
		}
	}
}

/**************************************************************************************************************************/
void DNAString::validate_char(char dna_char) const
{
	if((dna_char != 'A') && (dna_char != 'C') && (dna_char != 'G') && (dna_char != 'T') && (dna_char != 'N')) 
	{
		throw Illegal_DNAChar(dna_char, __PRETTY_FUNCTION__, "");
	}
}

/**************************************************************************************************************************/
DNAString DNAString::reverse_complement() const
{
	DNAString rc_str(::reverse(m_str));
	for(int i=rc_str.size()-1; i>=0; i--)
	{
		switch(rc_str.get(i)) 
		{
			case 'A':
				rc_str.set(i, 'T');
				break;
			case 'C':
				rc_str.set(i, 'G');
				break;
			case 'G':
				rc_str.set(i, 'C');
				break;
			case 'T':
				rc_str.set(i, 'A');
				break;
			// Else - N, leave as is
		}
	}
	return rc_str;
}

/***************************************************************************************************************************/
DNAString operator + (const DNAString& dna1, const DNAString& dna2)
{
        DNAString dna3 = dna1;

        dna3 += dna2;
        return dna3;
}

/***************************************************************************************************************************/
DNAString operator + (const DNAString& dna1, char ch)
{
        DNAString dna3 = dna1;

        dna3 += ch;
        return dna3;
}

/***************************************************************************************************************************/
ostream& operator << (ostream& os, const DNAString& dna)
{
        os << string(dna);
        return os;
}

/***************************************************************************************************************************/
double gc(const DNAString& dna_str)
{
	size_t gc = 0, Ns = 0;

	DNAString::const_iterator end=dna_str.end();
	for(DNAString::const_iterator it=dna_str.begin(); it!=end; it++) 
	{
		if(*it == 'N')
			Ns++;
		else if((*it == 'C') || (*it == 'G'))
			gc++;
	}
	if(dna_str.size() == Ns) {
		return 0;
	}

	return double(gc)/(dna_str.size()-Ns);
}

/***************************************************************************************************************************/
size_t Ns(const DNAString& dna_str)
{
	size_t Ns = 0;

	DNAString::const_iterator end=dna_str.end();
	for(DNAString::const_iterator it=dna_str.begin(); it!=end; it++) 
	{
		if(*it == 'N')
			Ns++;
	}

	return Ns;
}

}
