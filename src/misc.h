#ifndef MISC_H
#define MISC_H

#include <stdlib.h>
#include <time.h>
#include <vector>
#include <stdio.h>

using namespace std;

const char* get_time();

namespace MonteCarlo {

size_t Random();
void shuffle(vector<size_t>& values);
void shuffle(size_t* values, size_t n);

class Distribution {
public:
//        Distribution(const size_t* array, size_t n);
	Distribution() : values((size_t*)calloc(16384, sizeof(size_t))), n(0), values_size(16384)	{}
	~Distribution()						{free(values);}
	double  		cdf(size_t val) const;
	const Distribution&	operator += (size_t val);
	size_t			min_val() const			{return values[0];}
	size_t			max_val() const			{return values[n-1];}
	size_t  		size() const			{return n;}
protected:
        size_t*		values;
        size_t          n;
	size_t		values_size;
};

}

#endif // MISC_H
