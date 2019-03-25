#include "kdtree.h"

extern int samplingFactor;

int main(int argc, char **argv) {
	samplingFactor = atoi(argv[2]);
	parallelKDTree(argv[1], argv[3]);
	return 0;
}