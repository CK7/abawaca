#include "misc.h"
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include "TinyMT-src-1.0.1/tinymt/tinymt32.h"

#include<iostream>
using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const char* get_time()
{
	static int epoch = time(NULL);
	static char str[128];
	int t = time(NULL) - epoch;

	sprintf(str, "%02i%c%02i%c%02i", t/3600, ':', (t%3600)/60, ':', t%60);

	return str;
}

namespace MonteCarlo {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t Random()
{
        static tinymt32_t* tinymt = NULL;                               // Random Number Geberator
        if(tinymt == NULL) {
                tinymt = new tinymt32_t();
                tinymt->mat1 = strtoul("8f7011ee", NULL, 16);
                tinymt->mat2 = strtoul("fc78ff1f", NULL, 16);
                tinymt->tmat = strtoul("3793fdff", NULL, 16);
                srand(time(NULL));
                tinymt32_init(tinymt, rand());
        }
        return(tinymt32_generate_uint32(tinymt));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the Fisher-Yates shuffle algorithm
void shuffle(vector<size_t>& values)
{
	for(size_t i=values.size()-1; i>0; i--) {
		size_t j = Random()%i;
		size_t t = values[i];
		values[i] = values[j];
		values[j] = t;
	}  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the Fisher-Yates shuffle algorithm
void shuffle(size_t* values, size_t n)
{
	while(--n > 0) {
		size_t j = Random()%n;
		size_t t = values[n];
		values[n] = values[j];
		values[j] = t;
	}  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const Distribution& Distribution::operator += (size_t val)
{
	// Need to enlarge values
	if(n+1>values_size) {
cerr << "Reallocating" << endl;
		values = (size_t*)realloc(values, 2*values_size*sizeof(size_t));
		values_size *= 2;
cerr << "Done reallocating" << endl;
	}
	
	// values is sorted and we want to keep it that way
	size_t* p;
	for(p=values+n-1; (p >= values) && (val < *p); p--)
		*(p+1) = *p;

	// If we stopped somewhere in the middle or the end, or if we stopped at the beginning but the first 
	// value is smaller than val - insert to *(p+1)
	if(p!=values || val>*p)
		*(p+1) = val;
	else
		*p = val;
	n++;
	
	return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Distribution::cdf(size_t val) const
{
	// Easy cases
	if(val < *values)
		return 0;
	else if(val >= values[n-1])
		return 1;

	// We have to search distribution (modified binary search)
	size_t	imin=0, imax=n-1;
	while (imax >= imin)
	{
		/* calculate the midpoint for roughly equal partition */
      		int imid = (imin + imax) / 2;
 
      		// determine which subarray to search
      		if(values[imid] <= val)
		        imin = imid + 1;
     		 else if(values[imid] > val)
		        imax = imid - 1;
	}
	
	// imin must point to the variable before the one we are looking for
	return double(imin)/n;
}

}
