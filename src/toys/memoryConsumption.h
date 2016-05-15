#pragma once

#include <cstring>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atoi */

int memoryConsumptionParseLine(char* line){
        int i = std::strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
        i = atoi(line);
	return i;
}


int memoryConsumptionKB(){ //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = memoryConsumptionParseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

