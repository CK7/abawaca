/*
 * common.cpp
 *
 *  Created on: 31/Jan/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include "common.h"
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>

using namespace std;

time_t	epoch = time(NULL);

/****************************************************************************************************************/
string	time_since_epoch()
{
	char			str[32];
	unsigned int	t = time(NULL)-epoch;
	sprintf(str, "[%02d:%02d:%02d]\t", (unsigned int)(t/3600), (unsigned int)((t%3600)/60), (unsigned int)(t%60));
	return str;
}

/****************************************************************************************************************/
void chomp(string& line)
{
	string::reverse_iterator it = line.rbegin();
	while((it != line.rend()) && ((*it == '\n') || (*it == '\r'))) {
		line.erase(line.size()-1);
		it = line.rbegin();
	}
}

/****************************************************************************************************************/
void chomp(char* line)
{
	char* p = line+strlen(line);
	if(p == line)
		return;

	for(p--; ((*p == '\n') || (*p == '\r')); p--) {
		*p = 0;
		if(p == line)
			break;
	} 
}

/****************************************************************************************************************/
void remove_white_spaces(string& line)
{
	string::iterator it1 = line.begin(), it2 = line.begin();
	while(it2 != line.end()) {
		if(!isspace(*it2))
			*(it1++) = *it2;
		it2++;
	}
	line.erase(it2, line.end());
}

/****************************************************************************************************************/
void trim_front(string& line)
{
	string::iterator it1 = line.begin(), it2 = line.begin();
	if(it2 == line.end())
		return;
	if(!isspace(*(it2++)))
		return;

	while((it2 != line.end()) && isspace(*it2))
		it1++, it2++;

	line.erase(line.begin(), it1);
}

/****************************************************************************************************************/
void trim_back(string& line)
{
	string::reverse_iterator	it2 = line.rbegin();
	size_t			start_erase = line.size()-1;
	if(it2 == line.rend())
		return;
	if(!isspace(*(it2++)))
		return;

	while((it2 != line.rend()) && isspace(*it2))
		start_erase--, it2++;

	line.erase(start_erase);
}

/****************************************************************************************************************/
long split(char delimiter, string line, vector<string>& destination)
{
	destination.clear();
	if(line.size() == 0)
		return 0;

	size_t p1 = 0, p2 = 0;
	while((p1 < line.size()) && ((p2 = line.find(delimiter, p1)) != string::npos)) {
		destination.push_back(line.substr(p1, p2-p1));
		p1 = p2+1;
	}
	destination.push_back(line.substr(p1));
	return destination.size();
}

/****************************************************************************************************************/
string reverse(string str)
{
	string reverse_str;
	reverse_str.reserve(str.size()+1);
	for(string::reverse_iterator it = str.rbegin(); it != str.rend(); it++)
		reverse_str += *it;

	return reverse_str;
}

/****************************************************************************************************************/
int directory_exists(const char* dir_name)
{
	struct stat     buf;
	return((stat(dir_name, &buf) == 0) && S_ISDIR(buf.st_mode));
}

/****************************************************************************************************************/
int file_exists(const char* file_name)
{
        struct stat     buf;

        return((stat(file_name, &buf) == 0) && !S_ISDIR(buf.st_mode));
}

/****************************************************************************************************************/
bool read_line(string& line, FILE* fp)
{
	char str[1024];

	line.resize(0);
	while(fgets(str, 1024, fp)) {
		line += string(str);
		if(strchr(str, '\n'))
			break;
	}

	return !feof(fp);
}

/****************************************************************************************************************/
string reverse_complement(const string& dna_seq) 
{
	string rc_dna_seq = reverse(dna_seq);
	for(string::iterator it=rc_dna_seq.begin(); it!=rc_dna_seq.end(); it++) {
		switch(*it) {
			case 'A':
			case 'a':
				*it = 'T';
				break;
			case 'C':
			case 'c':
				*it = 'G';
				break;
			case 'G':
			case 'g':
				*it = 'C';
				break;
			case 'T':
			case 't':
				*it = 'A';
				break;
                }
        }
        return rc_dna_seq;
}

/****************************************************************************************************************/
uint64_t filesize(const char* filename)
{
	std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
	in.seekg(0, std::ifstream::end);
//	return (unsigned int64_t)(in.tellg()); 
	return (int64_t)(in.tellg()); 
}

/****************************************************************************************************************/
void print_percent(size_t old_p, size_t new_p, ostream& os)
{
	if(new_p == 0) {
		os << "0 %";
		return;
	}

	os << char(8) << char(8);
	if(old_p == 0)
		os << char(8);
	else {
		while(old_p > 0) {
			os << char(8);
			old_p /= 10;
		}
	}
	os << new_p << " %";
}

/****************************************************************************************************************/
vector<string> glob(const string& pat)
{
	glob_t glob_result;
	glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
	vector<string> ret;

	for(unsigned int i=0;i<glob_result.gl_pathc;++i){
		ret.push_back(string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);

	return ret;
}
