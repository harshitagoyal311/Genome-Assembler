#include "utils.h"

data *reservoirSampling2(dataArray *myArrayPointer, int k) {
	srand(time(NULL));
	data *tempDataList = NULL, *myDataList = myArrayPointer->dataList;
	int i, j, numCords = myArrayPointer->d, noOfPoints = myArrayPointer->n, kSelected[k];

	tempDataList = allocDataList(k, numCords);
	
	for(i = 0; i < k; i++)
		kSelected[i] = i;

	for(; i < noOfPoints; i++) {
		j = rand()%(i+1);
		if(j >= 0 && j < k)
			kSelected[j] = i;
	}

	for(i = 0; i < k; i++) {
		for(j = 0; j < numCords; j++) {
			tempDataList[i][j] = myDataList[kSelected[i]][j];
		}
	}

	return tempDataList;
}