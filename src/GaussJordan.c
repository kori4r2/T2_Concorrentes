/*
 *
// para compilar: mpicc GaussJordan.c -o GaussJordan -Wall -lm -fopenmp
// para rodar: mpirun -np 2 GaussJordan
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

#define NUM_THREADS 2
#define MAX(x,y) x>y?x:y

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
    double *pivo_row;
    struct pivo{
        double val;
        int ind;
    } local_pivo, pivo_reduce;

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
	ind_first_row = ((order/num_proc) )*my_rank +(my_rank<remainder?my_rank:remainder);
	row_status = (char *)calloc(sizeof(char), num_rows);
	pivo_row = (double *)malloc(sizeof(double) * order);
    
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
    printf("terminei scatter (rank %d)\n", my_rank);	

	int maxLoop = num_elements/order;
	for(int i = 0; i < maxLoop; i++){
		char string[256];
		sprintf(string, "(My rank is %d) | ", my_rank);
		for(int j = 0; j < order; j++){
			//printf("%s\n", string);
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

    //for(int i = 0; i < order; i++){
    for(int i = 0; i < 1; i++){
    	//verificar localmente qual o pivo -> struct (valor absoluto, idx)
    		//como verificar o indice?-
    		//como verificar qual indice(linha) ja foi usado?
	    		//->> numero de linhas = num_elements/order
	    		//->INDICE DA PRIMEIRA LINHA: ((order/num_proc) )*my_rank +(my_rank<remainder?my_rank:remainder)
	    		//VETOR com estado de cada linha
    	local_pivo.val = -1;
    	local_pivo.ind = -1;

    	for(int j = 0; j < num_rows; j++){

    		if(row_status[j] == 0){
    			//printf("(rank=%d) linha %d valida\n", my_rank, j);
    			if( fabs(recvbuf[j*order + pivo_col]) > local_pivo.val ){
    				local_pivo.val = fabs(recvbuf[j*order + pivo_col]);
    				local_pivo.ind = ind_first_row + j;
    			}
    		}

    	}

    	printf("eu rank = %d  | pivo local da col %d = %lf, linha %d\n", my_rank, pivo_col, local_pivo.val, local_pivo.ind);
    	//achar pivo global -> reduce (MAXLOC)

    	MPI_Reduce( &local_pivo, &pivo_reduce, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD ); 

    	if(my_rank == 0){
    		printf("(rank=%d) resultado do reduce: (%lf, %d)\n", my_rank, pivo_reduce.val, pivo_reduce.ind);
    	}

    	//broadcast rank do pivo -> marcar linha como usada

    	MPI_Bcast(&pivo_reduce, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);

    	//atualizar linha do pivo

    	int local_pivo_row = pivo_reduce.ind - ind_first_row;

    	int aux = pivo_reduce.ind / ((order/num_proc) + 1);
		int pivo_rank = aux < remainder? aux : remainder + ( (pivo_reduce.ind - ((order/num_proc) + 1) * remainder) / (order/num_proc));
    	

    	if(pivo_reduce.ind >= ind_first_row && pivo_reduce.ind <= ind_first_row+num_rows){
    		row_status[local_pivo_row] = 1;

    		for(int j = 0; j < order; j++){
	    		pivo_row[j] = recvbuf[local_pivo_row * order + j];
	    	}

	    	#pragma omp parallel num_threads(NUM_THREADS) shared (order)
			{
	    		//pragama parallel
	    			//pragam for
	    			//
	    		//
	    		int val = pivo_reduce.val;
	    		#pragma omp for shared(val, pivo_row, pivo_col)
	    		{
	    			for (int k = pivo_col+1; k < order; k++){
	    				pivo_row[k] = pivo_row[k]/pivo_row[pivo_col];
	    			}
	    		}
	    		pivo_row[pivo_col] = 1;
			}


			for(int j = 0; j < order; j++)
    			printf("(rank %d) linha pivo att: (%d) %lf\n", my_rank, j, pivo_row[j]);
    	}


    	//broadcast da linha

    	MPI_Bcast(pivo_row, order, MPI_DOUBLE, pivo_rank, MPI_COMM_WORLD);

    	printf("(rank=%d)bcast pivo row att: %lf\n", my_rank, pivo_row[2]);

    	//atualizar demais linhas

    		//pragama parallel
    			//pragam for
    			//
    		//

    	#pragma omp parallel num_threads(NUM_THREADS) shared (num_rows)
		{
    		//pragama parallel
    			//pragam for
    			//
    		//
    		for(int j = 0; j < num_rows; j++){
	    		if( ind_first_row + j != pivo_reduce.ind ){	
	    			double factor = recvbuf[j * order + pivo_col];
	    			printf("(rank %d) factor = %lf, row = %d\n",my_rank, factor, ind_first_row+j );
		    		#pragma omp for shared(recvbuf, pivo_row, pivo_col, factor, j, order)
		    		{
		    			for (int k = pivo_col+1; k < order; k++){
		    				recvbuf[j * order + k] = recvbuf[j * order + k] - (factor * pivo_row[k]);
		    				printf("===============(rank %d, thread %d)new value of (%d, %d): %lf\n", my_rank, omp_get_thread_num(), j+ind_first_row, k, recvbuf[j * order + k]);
		    			}
		    		}
		    		recvbuf[j * order + pivo_col] = 0;
		    	}
	    	}
		}



    	//atualizar pivo_col
    	pivo_col++;
    }

    	//VETOR COM A ORDENS DAS LINHAS DOS PIVOS
    //gather -> resultado

    //print resultado

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

