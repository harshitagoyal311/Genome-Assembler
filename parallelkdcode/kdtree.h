#include "utils.h"

int partition(dataArray *myArrayPointer, int low, int high, int dim);
double *quickselect(dataArray *myArrayPointer, int left, int right, int k, int dim);
void calcMinMaxOfDims(dataArray *myArrayPointer, int left, int right, int numCords, double *min, double *max);
int findMaxExtentDim(double *min, double *max, int numCords);
splitStruct *getSplitIndexArray(int nodes);
int calcMaxLevel(int nodes);
void printSplits(splitStruct *splitIndexArray, int splits);
void buildKDTree(dataArray *myArrayPointer, int left, int right, int currentLevel, int maxLevel, splitStruct *splitIndexArray, int *splitIndexPos);
void KDTree(char *in_file, int nodes, char *out_file);