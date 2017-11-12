/*
 *
// para compilar: mpicc GaussJordan.c -o GaussJordan -Wall -lm -fopenmp
// para rodar: mpirun -np 2 GaussJordan
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mpi.h"
#include <omp.h>

int getOrder(char * line, int length);
char * ReadLine(FILE *input, int * dim);
double * getMatrix(char * filename, int * dim );
double * getVector (char * filename, int size);
void destroyMatrix(double * vector, int order);

int main(int argc, char *argv[]){
    int *result;
    int order;
    double * matrix = NULL;
    double * vector=NULL;
    int my_rank, num_proc;      
    double * recvbuf;
    int * sendcounts = NULL, * displs = NULL;
    int pivo_col;

    int num_rows;
    int ind_first_row;
    char *row_status; //1-ja foi usada como pivo, 0-nao foi usada como pivo
    struct pivo{
        double val;
        int ind;
    } local_pivo;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

/*A-ler arquivos: matrix e vetor*/
    
    if(my_rank == 0){
        printf("sou rank 0\n");
        matrix = getMatrix("matriz.txt", &order);
		printf("terminei de ler a matriz\n");
		for(int i = 0; i < order; i++){
			for(int j = 0; j < order; j++){
				printf("%2.1lf  ", matrix[(i * order) + j]);
			}
			printf("\n");
		}
        vector = getVector("vetor.txt", order);
		printf("terminei de ler o vetor\n");
    }

	printf("Sou rank %d | antes:  order = %d\n", my_rank, order);
	MPI_Bcast(&order, 1, MPI_INT, 0, MPI_COMM_WORLD);
	printf("Sou rank %d | depois: order = %d\n", my_rank, order);
    	
/*B-scatter da matriz/vetor (varias linhas por processo)*/    
    
	printf("num_proc = %d\n", num_proc);
    int remainder = (order % num_proc);
    int num_elements = ((order / num_proc) + ((my_rank < remainder)? 1: 0)) * order;
	printf("Sour rank %d, num_elements = %d, remainder = %d, order = %d\n", my_rank, num_elements, remainder, order);

	num_rows = num_elements/order;
	ind_first_row = ((order/num_proc) )*my_rank +(my_rank<rest?my_rank:rest);
	row_status = (char *)calloc(sizeof(char), num_rows);
    
    recvbuf = (double *)malloc(sizeof(double) * (num_elements));
    
    if(my_rank == 0){
        sendcounts = (int *)malloc(sizeof(int) * num_proc);
        displs = (int *)malloc(sizeof(int) * num_proc);
            
        int sum = 0;
        
        for(int i = 0; i < num_proc; i++){
    		sendcounts[i] = ((order / num_proc) + ((i < remainder)? 1: 0)) * order;
			printf("sendcounts[%d] = %d\n", i, sendcounts[i]);
            
            displs[i] = sum;
            sum += sendcounts[i];
        }
    }
    
    MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, recvbuf, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	int maxLoop = num_elements/order;
	for(int i = 0; i < maxLoop; i++){
		char string[256];
		sprintf(string, "(My rank is %d) | ", my_rank);
		for(int j = 0; j < order; j++){
			sprintf(string, "%lf", recvbuf [(i * order) + j]);
		}
		printf("%s\n", string);
	}


/*C-Gauss jordan*/

    /*Passo 1: encontrar pivot */

    /*Passo 2: trocar linhas*/

    /*Passo 3: dividir elementos da linha k pelo pivot - n tarefas */

    /*Passo 4: somar linha k com os valores das demais linhas */
    
    pivo_col = 0;

    for(int i = 0; i < order; i++){
    	//verificar localmente qual o pivo -> struct (valor absoluto, idx)
    		//como verificar o indice?-
    		//como verificar qual indice(linha) ja foi usado?

	    		//->> numero de linhas = num_elements/order
	    		//->INDICE DA PRIMEIRA LINHA: ((order/num_proc) )*my_rank +(my_rank<rest?my_rank:rest)
	    		//VETOR com estado de cada linha
    	for(int j = 0; j < num_rows; j++){

    		

    	}


    	//achar pivo global -> reduce (MAXLOC)

    	//broadcast rank do pivo -> marcar linha como usada

    	//atualizar linha do pivo

    	//broadcast da linha

    	//atualizar demais linhas

    	//atualizar pivo_col
    	pivo_col++;
    }


    free(recvbuf);
    if(sendcounts != NULL)
    	free(sendcounts);
    if(displs != NULL)
    	free(displs);
    
    destroyMatrix(matrix, order);
    free(vector);

    MPI_Finalize();
    return 0;
}


void destroyMatrix(double * matrix, int order){
    free(matrix);
}

char * ReadLine(FILE * fp, int * dim){
    char * line = NULL;
    char c;
    int count =0;
    
    do {    
        c = fgetc(fp);
        
        if(c!='\n' && c!=EOF){
            line = (char*)realloc(line, sizeof(char)*(count+1));
            line [count++] = c;
             
        }else if(count>0){
            line = (char*)realloc(line,sizeof(char)*(count+1));
            line [count] = '\0';   
        
        }
    
    }while(c!=EOF && c!='\n');
 
    *dim = count;
    return line;
}


int getOrder(char * line, int length){
	int i, count = 0, readingNumber = 0;
	for(i = 0; i < length; i++){
		if(!readingNumber && isdigit(line[i])){
			readingNumber = 1;
			count++;
		}else if(readingNumber && !isdigit(line[i]) && line[i] != '.'){
			readingNumber = 0;
		}
	}
    return count;
}

double * getMatrix(char * filename, int * dim){

    char * line =NULL;

    double * matrix = NULL;

    int i ,k = 0,order = 0, length=0;

    FILE * fp = fopen(filename, "r");

    char * p, *q;

    
    
    line = ReadLine(fp, &length);
    
    order = getOrder (line, length);

    *dim = order;

    matrix = (double*)calloc(sizeof(double), order * order);
      
    
        
    do {
         
       p=&line[0]; 
      
       i=0;
       
       while(i<order){
    
            matrix[(k * order) + i]=strtod(p, &q);
			i++;
            q++;
        p = q;
                    
       }
       
       k++;
       
       free(line);
      
       line = ReadLine(fp, &length);
       
    
    }while(line!=NULL); 


    fclose(fp);
    return matrix;
}
double * getVector (char * filename, int size){
        FILE * fp = fopen(filename, "r");
		if(fp == NULL){
			printf("invalide filename\n");
			return NULL;
		}
        int i=0;
        double * vector = (double*)calloc(sizeof(double), size);

        while(i<size){
                fscanf(fp, "%lf", &vector[i]);
                i++;
        }    
        fclose(fp);
        return vector;
}

