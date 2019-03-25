#include "stdio.h"

void breakfile(char *filename, int noOfFiles, char *outfile) {
	FILE *fp = fopen(filename, "r"), *fp1 = NULL;
	int i, j, noOfPoints, numCords, pointsFile;
	fscanf(fp, " %d ", &noOfPoints);
	fscanf(fp, " %d ", &numCords);
	char outCycle[100], buffer[1000];
	pointsFile = noOfPoints/noOfFiles;
	int modPoint = noOfPoints % noOfFiles;
	//printf("Mod Point Is %d \n",modPoint);

	
	for(i = 0; i < noOfFiles; i++) {
		sprintf(outCycle, "%s%d", outfile, i);
		if(modPoint == 0){
			fp1 = fopen(outCycle, "w");
			fprintf(fp1, "%d\n%d\n", pointsFile, numCords);
			for(j = 0; j < pointsFile; j++) {
				fgets(buffer, 1000, fp);
				fputs(buffer, fp1);
			}
		}else{
			modPoint--;
			//printf("anand\n");
			int temp = pointsFile + 1;
			fp1 = fopen(outCycle, "w");
			fprintf(fp1, "%d\n%d\n", temp, numCords);
			for(j = 0; j < temp; j++) {
				fgets(buffer, 1000, fp);
				fputs(buffer, fp1);
			}
		}
		fclose(fp1);
	}
	

	fclose(fp);
}

int main(int argc, char **argv) {
	if(argc != 4) {
		printf("Arguments: i/p_file nodes o/p_file\n");
		return 0;
	}
	breakfile(argv[1], atoi(argv[2]), argv[3]);
}
