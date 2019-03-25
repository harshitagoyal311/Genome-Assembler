#include "kdtree.h"
#include <mpi.h>

dataArray *distributeDataSample(dataArray *myArrayPointer, samplingFunc sample, int sampleSize, int worldRank, int worldSize, MPI_Comm myComm);
int redistributeData(dataArray *myArrayPointer, int worldRank, int worldSize, int midValIndex, MPI_Comm myComm);
void buildParallelKDTreeTrueMedian(char *in_file);
void buildParallelKDTreeRandomSampling(dataArray *myArrayPointer, MPI_Comm myComm);
void parallelKDTree(char *in_file, char *out_file);