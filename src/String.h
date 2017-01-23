/*
 * String.h
 *
 *  Created on: 05/Feb/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef STRING_H_
#define STRING_H_

#include <iostream>

#include <string>
#include "common.h"
#include "Sequence.h"

using namespace std;

namespace Bio {

/***************************************************************************************************************************
 * Bio::DNAString
 * 
 ***************************************************************************************************************************/
class DNAString {
public:
	class Illegal_DNAString;
	class Illegal_DNAChar;
	typedef string::iterator iterator;
	typedef string::const_iterator const_iterator;
protected:
	string m_str;
public:
	DNAString() 									{}
	DNAString(const string& dna_string);
	const DNAString&	operator = (const string& s);
	size_t			size() const						{return m_str.size();}
	DNAString		subseq(size_t start = 0, size_t end = string::npos) const	{return (end != string::npos)? m_str.substr(start, (end-start+1)) : m_str.substr(start, string::npos);}
	const DNAString&	operator += (const DNAString& dna_str)			{m_str += dna_str.m_str; return *this;}
	const DNAString&	operator += (char dna_char) 				{this->validate_char(dna_char); m_str += toupper(dna_char); return *this;}
	operator string() const								{return m_str;}
	char			get(size_t i) const					{return m_str.at(i);}
	// Only const_iterator is allowed - we don't want someone messing up with m_str directly through iterator!
	const_iterator		begin() const						{return m_str.begin();}
	const_iterator		end() const						{return m_str.end();}
	iterator		begin()							{return m_str.begin();}
	iterator		end()							{return m_str.end();}
	const DNAString&	set(size_t i, char ch) 					{this->validate_char(ch); m_str.at(i) = ch; return *this;}
	DNAString		reverse_complement() const;
	bool			operator < (const DNAString& other) const		{return (m_str.compare(other.m_str) < 0);}
	bool			operator > (const DNAString& other) const		{return (other < *this);}
	bool			operator == (const DNAString& other) const		{return (!(*this < other) && !(other < *this));}
protected:
	void            	validate_string();
	void            	validate_char(char dna_char) const;
};

typedef Sequence<DNAString>	DNASequence;

DNAString 	operator + (const DNAString& dna1, const DNAString& dna2);
DNAString 	operator + (const DNAString& dna1, char ch);
ostream& 	operator << (ostream& os, const DNAString& dna);
double		gc(const DNAString& dna_str);
size_t		Ns(const DNAString& dna_str);

/***************************************************************************************************************************/
class DNAString::Illegal_DNAString : public exception_base {
public:
        Illegal_DNAString(string str, const char* func, string err_msg) : exception_base("", "")
                {
                        msg = string("Fatal error, attempted to initialize DNAString with illegal string:\n\n") + str;
                        if(err_msg.size() != 0)
                                msg += string("\n\nError message:\t") + err_msg + string("\n");
			else
				msg += "\n\n";
                        msg += string("Function:\t") + string(func) + string("\n");
                }
};

/***************************************************************************************************************************/
class DNAString::Illegal_DNAChar : public exception_base {
public:
        Illegal_DNAChar(char ch, const char* func, string err_msg) : exception_base("", "")
                {
			char str[2] = {ch, 0};
                        msg = string("Fatal error, attempted to initialize DNAString with illegal character:\n\n") + string(str);
                        if(err_msg.size() != 0)
                                msg += string("\n\nError message:\t") + err_msg + string("\n");
			else
				msg += "\n\n";
                        msg += string("Function:\t") + string(func) + string("\n");
                }
};

}

#endif /* STRING_H_ */
