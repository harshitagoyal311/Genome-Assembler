#include "parallelkdtree.h"
#include <mpi.h>

int samplingFactor;

dataArray *distributeDataSample(dataArray *myArrayPointer, samplingFunc sample, int sampleSize, int worldRank, int worldSize, MPI_Comm myComm) {
	int i, j, dataSize = 0, numCords = myArrayPointer->d, sampleSizeArray[worldSize], displs[worldSize];
	data *sendData = sample(myArrayPointer, sampleSize);

	MPI_Allgather(&sampleSize, 1, MPI_INT, sampleSizeArray, 1, MPI_INT, myComm);

	displs[0] = 0;
	sampleSizeArray[0] = sampleSizeArray[0]*numCords;
	for(i = 1; i < worldSize; i++) {
		sampleSizeArray[i] = sampleSizeArray[i]*numCords;
		displs[i] = displs[i-1]+sampleSizeArray[i-1];
	}

	for(i = 0; i < worldSize; i++)
		dataSize += sampleSizeArray[i];
	dataSize = dataSize/numCords;

	data *recvData = allocDataList(dataSize, numCords);

	MPI_Allgatherv(sendData[0], sampleSize*numCords, MPI_DOUBLE, recvData[0], sampleSizeArray, displs, MPI_DOUBLE, myComm);
	MPI_Barrier(myComm);

	freeDataList(sendData);
	dataArray *recvDataPointer = (dataArray *)malloc(sizeof(dataArray));
	recvDataPointer->n = dataSize;
	recvDataPointer->d = numCords;
	recvDataPointer->dataList = recvData;
	return recvDataPointer;
}

int redistributeData(dataArray *myArrayPointer, int worldRank, int worldSize, int midValIndex, MPI_Comm myComm) {
	data *myDataList = myArrayPointer->dataList;
	int numCords = myArrayPointer->d, noOfPoints = myArrayPointer->n, source = (worldRank+worldSize/2)%worldSize, destination, recvCount, sendCount, newCount;
	destination = source;
	MPI_Status recvStatus;

	if(worldRank >= worldSize/2) {
		sendCount = midValIndex+1;
		MPI_Probe(source, 0, myComm, &recvStatus);
		MPI_Get_count(&recvStatus, MPI_DOUBLE, &recvCount);
		recvCount = recvCount/numCords;
		newCount = noOfPoints - sendCount + recvCount;
		data *newDataList = allocDataList(newCount, numCords);
		memcpy(newDataList[0], myDataList[sendCount], sizeof(double)*numCords*(noOfPoints-sendCount));
		MPI_Recv(newDataList[noOfPoints-sendCount], recvCount*numCords, MPI_DOUBLE, source, 0, myComm, MPI_STATUS_IGNORE);
		//MPI_Sendrecv(myDataList[0], sendCount*numCords, MPI_DOUBLE, destination, 0, newDataList[noOfPoints-sendCount], recvCount*numCords, MPI_DOUBLE, source, 0, myComm, MPI_STATUS_IGNORE);
		MPI_Send(myDataList[0], sendCount*numCords, MPI_DOUBLE, destination, 0, myComm);
		myArrayPointer->dataList = newDataList;
		freeDataList(myDataList);
		myArrayPointer->n = newCount;
		//printf("process: %d, newcount = %d\n", worldRank, newCount);
		return newCount;
	}

	else {
		sendCount = noOfPoints - (midValIndex + 1);
		MPI_Send(myDataList[midValIndex+1], sendCount*numCords, MPI_DOUBLE, destination, 0, myComm);
		MPI_Probe(source, 0, myComm, &recvStatus);
		MPI_Get_count(&recvStatus, MPI_DOUBLE, &recvCount);
		recvCount = recvCount/numCords;
		newCount = noOfPoints - sendCount + recvCount;
		data *newDataList = allocDataList(newCount, numCords);
		memcpy(newDataList[0], myDataList[0], sizeof(double)*numCords*(midValIndex+1));
		MPI_Recv(newDataList[midValIndex+1], recvCount*numCords, MPI_DOUBLE, source, 0, myComm, MPI_STATUS_IGNORE);
		//MPI_Sendrecv(myDataList[midValIndex+1], sendCount*numCords, MPI_DOUBLE, destination, 0, newDataList[midValIndex+1], recvCount*numCords, MPI_DOUBLE, source, 0, myComm, MPI_STATUS_IGNORE);
		myArrayPointer->dataList = newDataList;
		freeDataList(myDataList);
		myArrayPointer->n = newCount;
		//printf("process: %d, newcount = %d\n", worldRank, newCount);
		return newCount;
	}

	return;
}
/*
void buildParallelKDTreeTrueMedian(char *in_file) {
	MPI_Init(NULL, NULL);

	int worldSize, worldRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	char s[100];
	sprintf(s, "%s%d", in_file, worldRank);
	int numCords, noOfPoints, i, maxLevel, maxExtentDim;
	maxLevel = calcMaxLevel(worldSize);
	dataArray myArray = readDataFromFile(s);
	dataArray *myArrayPointer = &myArray;
	numCords = myArrayPointer->d;
	noOfPoints = myArrayPointer->n;

	double min[numCords], max[numCords];

	calcMinMaxOfDims(myArrayPointer, 0, noOfPoints-1, numCords, min, max);

	//printf("Local min-max for process %d:\n", worldRank);
	// for(i = 0; i < numCords; i++)
	// 	printf("dim: %d, min: %lf, max: %lf\n", i, min[i], max[i]);

	MPI_Allreduce(MPI_IN_PLACE, min, numCords, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, max, numCords, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);	

	MPI_Barrier(MPI_COMM_WORLD);

	maxExtentDim = findMaxExtentDim(min, max, numCords);
	// printf("Global min-max for process %d:\n", worldRank);
	// for(i = 0; i < numCords; i++)
	// 	printf("dim: %d, min: %lf, max: %lf\n", i, min[i], max[i]);

	MPI_Finalize();
}
*/

void buildParallelKDTreeRandomSampling(dataArray *myArrayPointer, MPI_Comm myComm) {
	int worldSize, worldRank;
	MPI_Comm_rank(myComm, &worldRank);
	MPI_Comm_size(myComm, &worldSize);

	if(worldSize == 1) {
		MPI_Finalize();
		return;
	}

	int color, numCords, noOfPoints, i, maxLevel, maxExtentDim, midValIndex, sampleMidValIndex, sampleSize;
	data *myDataList = myArrayPointer->dataList;
	numCords = myArrayPointer->d;
	noOfPoints = myArrayPointer->n;

	samplingFunc mySampleFunc = &reservoirSampling2;
	sampleSize = noOfPoints/samplingFactor;
	dataArray *sampleArrayPointer = distributeDataSample(myArrayPointer, mySampleFunc, sampleSize, worldRank, worldSize, myComm);
	// dataArray sampleArray;
	// sampleSize = sampleSize*worldSize;
	// sampleArray.n = sampleSize;
	// sampleArray.d = numCords;
	// sampleArray.dataList = sampleDataList;
	// dataArray *sampleArrayPointer = &sampleArray;
	sampleSize = sampleArrayPointer->n;
	data *sampleDataList = sampleArrayPointer->dataList;
	double min[numCords], max[numCords], sampleMidVal;
	calcMinMaxOfDims(sampleArrayPointer, 0, sampleSize-1, numCords, min, max);
	maxExtentDim = findMaxExtentDim(min, max, numCords);
	quickselect(sampleArrayPointer, 0, sampleSize-1, sampleSize/2+1, maxExtentDim);

	MPI_Barrier(MPI_COMM_WORLD);

	// for(i = 0; i < worldSize; i++) {
	// 	if(worldRank == i) {
	// 		printf("Process: %d\nSample:\n", worldRank);
	// 		printDataList(sampleDataList, sampleSize, numCords);
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }	
	
	sampleMidValIndex = sampleSize/2;
	printf("[+]Sample Size for Process %d/%d is %d\n", worldRank, worldSize, sampleSize);
	sampleMidVal = sampleDataList[sampleMidValIndex][maxExtentDim];
	MPI_Barrier(MPI_COMM_WORLD);
	if(worldRank == 0)
		printf("\n");
	//printf("%d: Midval = %lf, Dim = %d\n", worldRank, sampleMidVal, maxExtentDim);
	MPI_Barrier(MPI_COMM_WORLD);

	midValIndex = partitionAroundValue(myArrayPointer, 0, noOfPoints-1, maxExtentDim, sampleMidVal);

	// for(i = 0; i < worldSize; i++) {
	// 	if(worldRank == i) {
	// 		printf("Process: %d\nMedian Index: %d, Value: %lf, Dim: %d\nMain:\n", worldRank, midValIndex, sampleMidVal, maxExtentDim);
	// 		printDataList(myArray.dataList, noOfPoints, numCords);
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }

	redistributeData(myArrayPointer, worldRank, worldSize, midValIndex, myComm);
	MPI_Barrier(MPI_COMM_WORLD);

	for(i = 0; i < worldSize; i++) {
		if(worldRank == i) {
			//printf("Process: %d\nSize: %d, Value: %lf, Dim: %d", worldRank, myArrayPointer->n, sampleMidVal, maxExtentDim);
			//printDataList(myArrayPointer->dataList, myArrayPointer->n, numCords);
			//printf("Array Dim: %d, size: %d\n", myArrayPointer->d, myArrayPointer->n);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	color = (worldRank < worldSize/2) ? 0 : 1;
	MPI_Comm splitComm;
	MPI_Comm_split(myComm, color, worldRank, &splitComm);
	int newRank, newSize;
	MPI_Comm_rank(splitComm, &newRank);
	MPI_Comm_size(splitComm, &newSize);
	//printf("Proc: %d/%d to %d/%d\n", worldRank, worldSize, newRank, newSize);
	MPI_Barrier(myComm);
	buildParallelKDTreeRandomSampling(myArrayPointer, splitComm);
	return;
}

void parallelKDTree(char *in_file, char *out_file) {
	MPI_Init(NULL, NULL);
	int numCords, noOfPoints, worldRank;
	MPI_Comm myComm = MPI_COMM_WORLD;
	MPI_Comm_rank(myComm, &worldRank);
	char s[100];
	sprintf(s, "%s%d", in_file, worldRank);
	dataArray myArray = readDataFromFile(s);
	dataArray *myArrayPointer = &myArray;
	numCords = myArrayPointer->d;
	noOfPoints = myArrayPointer->n;

	buildParallelKDTreeRandomSampling(myArrayPointer, MPI_COMM_WORLD);
	sprintf(s, "%s%d", out_file, worldRank);
	writeDataToFile(*myArrayPointer, s);
	return;
}