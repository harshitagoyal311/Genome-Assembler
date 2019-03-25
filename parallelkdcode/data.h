#include "def.h"

extern clock_t START_TIME, END_TIME;
extern double CPU_TIME;

data *allocDataList(int noOfPoints, int numCords);
dataArray readDataFromFile(char *filename);
void printPoint(data myPoint, int numCords);
void printDataArray(dataArray myArray);
void printDataList(data *myDataList, int noOfPoints, int numCords);
void freeData(dataArray myArray);
void writeDataToFile(dataArray myArray, char *filename);
int findIndex(data *myDataList, double val, int valDim, int noOfPoints);
int power_2(int e);
void start_time();
void end_time();
void breakfile(char *filename, int noOfFiles, char *outfile);