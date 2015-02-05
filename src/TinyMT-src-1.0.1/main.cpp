#include "MonteCarlo.h"

int main(int argc, char * argv[]) 
{
	if(argc < 2) {
		fprintf(stderr, "Usage: %s <scaf dp dist-file>\n\n", argv[0]);
		exit(-1);
	}
	MonteCarlo	rng(argv[1]);
}
