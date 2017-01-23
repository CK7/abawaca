/*
 * common.h
 *
 *  Created on: 31/Jan/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <bitset>
#include <ctype.h>
#include <time.h>
#include <glob.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <inttypes.h>

using namespace std;

// Will return a string containing the time passed from the moment the program started running
string	time_since_epoch();

// Similar to perl's chomp - will trim any \n and \r from the end of the string
void chomp(string& line);
void chomp(char* line);

// Will remove any white space that is recognized by isspace() from the string
void remove_white_spaces(string& line);

// Will remove any white space that is recognized by isspace() from the beginning of the string until the first non-space char
void trim_front(string& line);

// Will remove any white space that is recognized by isspace() from the back of the string until the first non-space char
void trim_back(string& line);

// Simulates perl's split, will store result in desctination and return the number of fields read
long split(char separator, string line, vector<string>& destination);

// Reverse a string
string reverse(string line);

// Returns true if dir_name is an existing directory, false otherwise
int directory_exists(const char* dir_name);

// Returns true if file_name is an existing file, false otherwise
int file_exists(const char* file_name);

// Will return all files that fit the glob path provided by pat
vector<string> glob(const string& pat);

// Reads a line from a file, everything up to the newline sign. Returns true if can keep reading from fp, false if reached EOF
bool read_line(string& line, FILE* fp);

// Returns the reverse complement of the input DNA sequence
string reverse_complement(const string& dna_seq);

// Returns the size of a file
uint64_t filesize(const char* filename);

// Prints "new_p %" (e.g. "75 %" if new_p=75) in the same position on the screen where old_p used to be.
// This function assumes that the position on the screen did not change since old_p was written.
void print_percent(size_t old_p, size_t new_p, ostream& os=cerr);

#endif /* COMMON_H_ */
