#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

typedef struct dataArray dataArray;
typedef struct splitStruct splitStruct;
typedef double *data;

struct dataArray {
	data *dataList;
	int n;
	int d;
};

struct splitStruct {
	int val;
	int dim;
};