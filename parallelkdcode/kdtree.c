#include "kdtree.h"

int partition(dataArray *myArrayPointer, int left, int right, int dim) {
	int i = left, j, k, numCords = myArrayPointer->d;
	data *myDataList = myArrayPointer->dataList;
	double x = myDataList[right][dim], temp;
	for(j = left; j < right; j++) {
		if(myDataList[j][dim] <= x) {
			for(k = 0; k < numCords; k++) {
				temp = myDataList[j][k];
				myDataList[j][k] = myDataList[i][k];
				myDataList[i][k] = temp;	
			}
			i++;
		}
	}

	for(k = 0; k < numCords; k++) {
		temp = myDataList[right][k];
		myDataList[right][k] = myDataList[i][k];
		myDataList[i][k] = temp;	
	}
	return i;
}

int randomPartition(dataArray *myArrayPointer, int left, int right, int dim) {
	srand(time(NULL));
	int p = left + rand()%(right-left+1), numCords = myArrayPointer->d, k;
	data *myDataList = myArrayPointer->dataList;
	double temp;
	for(k = 0; k < numCords; k++) {
		temp = myDataList[right][k];
		myDataList[right][k] = myDataList[p][k];
		myDataList[p][k] = temp;	
	}
	usleep(100);

	return partition(myArrayPointer, left, right, dim);
}

int partitionAroundValue(dataArray *myArrayPointer, int left, int right, int dim, double val) {
	int i = left, j, k, numCords = myArrayPointer->d;
	data *myDataList = myArrayPointer->dataList;
	double temp;
	for(j = left; j <= right; j++) {
		if(myDataList[j][dim] <= val) {
			for(k = 0; k < numCords; k++) {
				temp = myDataList[j][k];
				myDataList[j][k] = myDataList[i][k];
				myDataList[i][k] = temp;	
			}

			i++;
		}
	}

	return i-1;
}

double *quickselect(dataArray *myArrayPointer, int left, int right, int k, int dim) {
	
	int pivotIndex = 0;

	while(k > 0 && k <= right-left+1) {
		pivotIndex = randomPartition(myArrayPointer, left, right, dim);

		if(pivotIndex-left == k-1)
			return myArrayPointer->dataList[pivotIndex];

		else if(pivotIndex-left > k-1) {
			right = pivotIndex-1;
			//return quickselect(myArrayPointer, left, pivotIndex-1, k, dim);
		}

		else {
			k = k-pivotIndex+left-1;
			left = pivotIndex+1;
			//return quickselect(myArrayPointer, pivotIndex+1, right, k-pivotIndex+left-1, dim);
		}
	}

	return NULL;
}

void calcMinMaxOfDims(dataArray *myArrayPointer, int left, int right, int numCords, double *min, double *max) {
	data *myDataList = myArrayPointer->dataList;
	int n = myArrayPointer->n;
	int i = 0, j = 0;
	double temp;

	for(j = 0; j < numCords; j++) {
		min[j] = myDataList[left][j];
		max[j] = min[j];
	}

	for(i = left; i <= right; i++) {
		for(j = 0; j < numCords; j++) {
			temp = myDataList[i][j];
			if(temp < min[j])
				min[j] = temp;
			if(temp > max[j])
				max[j] = temp;
		}
	}
}

int findMaxExtentDim(double *min, double *max, int numCords) {
	int i = 0, maxExtentIndex = 0;
	double maxExtent = 0, temp = 0;
	for(i = 0; i < numCords; i++) {
		temp = max[i]-min[i];
		if(temp > maxExtent) {
			maxExtent = temp;
			maxExtentIndex = i;
		}
	}

	return maxExtentIndex;
}

splitStruct *getSplitIndexArray(int nodes) {
	splitStruct *splitIndexArray = (splitStruct *)malloc(sizeof(splitStruct)*(nodes+1));
	splitIndexArray[0].dim = -1;
	splitIndexArray[0].val = 0;
	return splitIndexArray;
}

int calcMaxLevel(int nodes) {
	return ceil(log2(nodes));
}

void printSplits(splitStruct *splitIndexArray, int splits) {
	int i = 0;
	for(i = 0; i < splits; i++)
		printf("\nDIM: %d, SPLIT: %d", splitIndexArray[i].dim, splitIndexArray[i].val);
	printf("\n\n");
}

void buildKDTree(dataArray *myArrayPointer, int left, int right, int currentLevel, int maxLevel, splitStruct *splitIndexArray, int *splitIndexPos) {
	if(currentLevel == maxLevel)
		return;

	int numCords, splitIndex, maxExtentDim, noOfPoints, midValIndex;
	numCords = myArrayPointer->d;
	noOfPoints = myArrayPointer->n;
	splitIndex = (right-left+1)/2;
	double min[numCords], max[numCords];

	calcMinMaxOfDims(myArrayPointer, left, right, numCords, min, max);
	maxExtentDim = findMaxExtentDim(min, max, numCords);
	quickselect(myArrayPointer, left, right, splitIndex, maxExtentDim);
	midValIndex = left+splitIndex-1;
	printf("level: %d, dim: %d, midValIndex: %d\n", currentLevel, maxExtentDim, midValIndex);
	
	buildKDTree(myArrayPointer, left, midValIndex-1, currentLevel+1, maxLevel, splitIndexArray, splitIndexPos);
	
	splitIndexArray[*splitIndexPos].val = midValIndex;
	splitIndexArray[*splitIndexPos].dim = maxExtentDim;
	*splitIndexPos = *splitIndexPos + 1;
	
	buildKDTree(myArrayPointer, midValIndex+1, right, currentLevel+1, maxLevel, splitIndexArray, splitIndexPos);
	
	return;
}

void KDTree(char *in_file, int nodes, char *out_file) {
	int maxLevel = calcMaxLevel(nodes), splitIndexPos = 1;
	dataArray myArray = readDataFromFile(in_file);
	dataArray *myArrayPointer = &myArray;
	
	splitStruct *splitIndexArray = getSplitIndexArray(nodes);
	start_time();
	buildKDTree(myArrayPointer, 0, myArray.n-1, 0, maxLevel, splitIndexArray, &splitIndexPos);
	end_time();
	splitIndexArray[splitIndexPos].val = myArray.n;
	splitIndexArray[splitIndexPos].dim = -1;
	
	printSplits(splitIndexArray, nodes+1);
	assignNodes(myArray, splitIndexArray, nodes, out_file);
}