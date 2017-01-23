/*
 * exceptions.h
 *
 *  Created on: 31/Jan/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <stdio.h>
#include <exception>
#include "Sequence.h"

using namespace std;

//namespace Bio {

/****************************************************************************************************************/
class exception_base : public exception {
public:
	exception_base(string err_msg, const char* func) : msg(err_msg)
		{
			msg += string("\nFunction:\t") + string(func);
			msg += string("Contact itaish@berkeley.edu for more information");
		}
	~exception_base() throw()							{}
	virtual const char* what() const throw()			{return msg.c_str();}
protected:
	string msg;
};

namespace Bio {

/****************************************************************************************************************/
class Bad_file : public exception_base {
public:
	Bad_file(string illegal_file, string error_type, const char* func, size_t line_num = 0) : exception_base("", "")
		{
			char line[32];
			if(line_num) {
				sprintf(line, "%u", (unsigned int)line_num);
				msg = string("Fatal error, line ") + string(line) + string(" file ") + illegal_file + string(": ") + error_type;
			}
			else {
				msg = string("Fatal error, file ") + illegal_file + string(": ") + error_type;
			}
			msg += string("\nFunction:\t") + string(func) + string("\n");
		}
};

/****************************************************************************************************************/
class Missing_sequence : public exception_base {
public:
	Missing_sequence(string display_id, const char* func, string err_msg, string db_file) :
		exception_base(string("Fatal error, missing sequence:\ndisplay_id:\t")+display_id+string("\nDatabase:\t")+db_file+string("Error msg:\t")+err_msg, func)
		{}
};

}

/****************************************************************************************************************/
class Software_bug : public exception {
	string msg;
public:
	Software_bug(string description, string additional_info, const char* func /*=__PRETTY_FUNCTION__*/)
		{
			msg = string("Software bug, please report to itaish@berkeley.edu:\nDescription:\t") + description + string("\n");
			if(additional_info.size() > 0)
				msg += string("Information:\t") + additional_info + string("\n");
			if(func)
				msg += string("Bug occurred in function ") + string(func) + string("\n");
		}
	~Software_bug() throw()
		{}
	virtual const char* what() const throw()
		{return msg.c_str();}
};

/****************************************************************************************************************/
template <class S>
class Illegal_sequence : public exception_base {
public:
	Illegal_sequence(const S& seq, const char* func, string err_msg) : exception_base("", "")
	{
		msg = string("Fatal error, illegal sequence found:\ndisplay_id:\t") + seq.display_id() + string("\nDescription:\t") +
			seq.desc();
		if(err_msg.size() != 0)
		{
			msg += string("Error message:\t") + err_msg + string("\n");
		}
		msg += string("Function:\t") + string(func) + string("\n");
		
	}
};
//}
#endif /* EXCEPTIONS_H_ */
