#include "data.h"

clock_t START_TIME, END_TIME;
double CPU_TIME;

data *allocDataList(int noOfPoints, int numCords) {
	int i = 0;
    double *temp = (double *)malloc(noOfPoints*numCords*sizeof(double));
    double **tempDataList = (double **)malloc(noOfPoints*sizeof(double*));
    for(i=0; i<noOfPoints; i++)
        tempDataList[i] = &(temp[numCords*i]);

    return tempDataList;
}

dataArray readDataFromFile(char *filename) {
	FILE *fp = fopen(filename, "r");
	int i = 0, n = 0, d = 0, j = 0;

	dataArray myArray;
	fscanf(fp, " %d ", &myArray.n);
	fscanf(fp, " %d ", &myArray.d);

	n = myArray.n;
	d = myArray.d;
	data *tempDataList = allocDataList(n, d);

	for(i = 0; i < n; i++) {
		for(j = 0; j < d; j++) {
			fscanf(fp, " %lf", &tempDataList[i][j]);
		}
	}

	myArray.dataList = tempDataList;

	fclose(fp);
	return myArray;
}

void printPoint(data myPoint, int numCords) {
	int i = 0;
	for(i = 0; i < numCords; i++)
		printf("%lf ", myPoint[i]);
	printf("\n");
}

void printDataArray(dataArray myArray) {
	int i = 0, n = 0, d = 0, j = 0;
	n = myArray.n;
	d = myArray.d;

	for(i = 0; i < n; i++) {
		printf("%d: ", i+1);
		for(j = 0; j < d; j++) {
			printf("%lf ", myArray.dataList[i][j]);
		}

		printf("\n");
	}

	printf("\n");
}

void printDataList(data *myDataList, int noOfPoints, int numCords) {
	int i = 0, j = 0;

	for(i = 0; i < noOfPoints; i++) {
		printf("%d: ", i+1);
		for(j = 0; j < numCords; j++) {
			printf("%lf ", myDataList[i][j]);
		}

		printf("\n");
	}

	printf("\n");
}

void freeDataList(data *myDataList) {
	free(myDataList[0]);
	free(myDataList);
}

void writeDataToFile(dataArray myArray, char *filename) {
	FILE *fp = fopen(filename, "w");
	int i = 0, n = 0, d = 0, j = 0;
	n = myArray.n;
	d = myArray.d;
	fprintf(fp, "%d\n", n);
	fprintf(fp, "%d\n", d);
	for(i = 0; i < n; i++) {
		//fprintf(fp, "%d\t", i);
		for(j = 0; j < d; j++)
			fprintf(fp, "%lf ", myArray.dataList[i][j]);

		fprintf(fp, "\n");
	}
}

int power_2(int e) {
	int p = 1;
	while(e-- > 0) {
		p = p << 1;
	}
	return p;
}

void assignNodes(dataArray myArray, splitStruct *splitIndexArray, int nodes, char *filename) {
	int i = 0, j = 0, k = 0, d = myArray.d, n = myArray.n;
	FILE *fp = fopen(filename, "w");
	fprintf(fp, "%d\n%d\n%d\n", n, d, nodes);
	for(i = 1; i <= nodes; i++) {
		for(j = splitIndexArray[i-1].val; j < splitIndexArray[i].val; j++) {
			for(k = 0; k < d; k++)
				fprintf(fp, "%lf ", myArray.dataList[j][k]);
			fprintf(fp, "%d\n", i-1);
		}
	}
}

int findIndex(data *myDataList, double val, int valDim, int noOfPoints) {
	int i = 0; 
	for(i = 0; i < noOfPoints; i++) {
		if(myDataList[i][valDim] == val)
			return i;
	}
	return -1;
} 

void start_time() {
	START_TIME = clock();
}

void end_time() {
	END_TIME = clock();
	CPU_TIME = ((double) (END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("CPU TIME: %lf\n", CPU_TIME);
}