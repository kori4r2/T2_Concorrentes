#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

char *myReadLine(FILE *input, bool *endFileReached){
	char inChar, output = NULL;
	int count = 0;

	inChar = fgetchar(input);
	while(inChar != '\n' && inChar != EOF){
		output = (char*)realloc(output, sizeof(char) * (count + 1));
		output[count++] = inChar;
		inChar = fgetchar(input);
	}
	if(count == 0){
		fprintf(stderr, "Function called at EOF");
		return NULL;
	} else{
		output[count++] = '\0';
		if(inChar == EOF)
			(*endFileReached) = true;
	}
	return output;
}

int main(int argc, char *argv[]){
	int *result;

	return 0;
}
