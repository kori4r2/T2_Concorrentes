/*
 * Nome 								- NUSP
 * Henrique Caetano Anraki 				– 8643412
 * Italo Tobler Silva 					– 8551910
 * Marcus Vinicius dos Santos Araujo 	– 9005871
 *
 * Instrucoes para execucao...
 * Em uma maquina local, sem makefile:
 * 	Para compilar: mpicc GaussJordan.c -o GaussJordan -lm -fopenmp
 * 	Para rodar: mpirun -np 2 GaussJordan
 * Em uma maquina local, usando o makefile na raiz do projeto:
 * 	Para compilar: make
 * 	Para rodar: make mpi_run
 * 	Para enviar os arquivos necessarios ao cluster disponibilizado: make sendfiles
 * 	Para se conectar ao mesmo cluster usando ssh: make remote
 * No cluster remoto disponibilizado para a disciplina, usando o makefile na pasta remoteFiles:
 * 	Para compilar: make
 * 	Para rodar: make run
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

#define NUM_THREADS 8

/* A partir de uma linha do arquivo de entrada, recupera a ordem da matriz pela contagem de elementos*/
int getOrder(char * line, int length); 
/* Le um arquivo linha por linha*/
char * ReadLine(FILE *input, int * dim);
/* Recupera uma matriz e sua ordem a partir de um arquivo. Arquivo : .txt*/
double * getMatrix(char * filename, int * dim );
/* Dado um arquivo de entrada e a quantidade de linhas a serem computadas, retorna um vetor de elementos*/
double * getVector (char * filename, int size);

int main(int argc, char *argv[]){
    
    int i, j, k;
    int order;
    double *matrix = NULL;
    double *vector=NULL;
    int my_rank, num_proc;      
    double *recvbuf, *recvbuf_v;
    int *sendcounts = NULL, *displs = NULL;
    int pivo_col;

    int num_rows;
    int ind_first_row;
    char *row_status; //1-ja foi usada como pivo, 0-nao foi usada como pivo
    double *pivo_row;
    double pivo_row_v;
    struct pivo{
        double val;
        int ind;
    } local_pivo, pivo_reduce;

    int *pivo_order = NULL;
    double time_i = 0, time_f;


/*Inicia processos e comunicadores*/

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

/*Processo Pai : Le arquivos de entrada*/
    
    if(my_rank == 0){   
		matrix = getMatrix("matriz.txt", &order);
		vector = getVector("vetor.txt", order);
		pivo_order = (int *)malloc(sizeof(int) * order);
    }

	// Repassa o tamanho da matriz a todos os processos
    MPI_Bcast(&order, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	
/*Distribui (Scatter) linhas das matrizes pelos processos*/    

    /*calculo de carga ou granularidade e alocacao de buffer/indice*/ 
    int remainder = (order % num_proc); 
    int num_elements = ((order / num_proc) + ((my_rank < remainder)? 1: 0)) * order;

    num_rows = num_elements/order;
    ind_first_row = ((order/num_proc) )*my_rank +(my_rank<remainder?my_rank:remainder);
    row_status = (char *)calloc(sizeof(char), num_rows);
    pivo_row = (double *)malloc(sizeof(double) * order);
    
    recvbuf = (double *)malloc(sizeof(double) * (num_elements));
    recvbuf_v = (double *)malloc(sizeof(double) * (num_rows));
    
    if(my_rank == 0){
        sendcounts = (int *)malloc(sizeof(int) * num_proc);
        displs = (int *)malloc(sizeof(int) * num_proc);
            
        int sum = 0;
        
        for( i = 0; i < num_proc; i++){
    		sendcounts[i] = ((order / num_proc) + ((i < remainder)? 1: 0)) * order;
            displs[i] = sum;
            sum += sendcounts[i];
        }
    }

    /*Medicao do tempo de execucao*/
    if(my_rank == 0)
    	time_i = omp_get_wtime();
    

	/*Distribuicao dos valores pelos processos*/
    MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, recvbuf, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(my_rank == 0){
	    for(i = 0; i < num_proc; i++){
	    	sendcounts[i] /= order;
	    	displs[i] /= order;
	    }
	}

    MPI_Scatterv(vector, sendcounts, displs, MPI_DOUBLE, recvbuf_v, num_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	

/*Algoritmo de Gauss Jordan*/
    
	// Percorre todas as colunas da matriz em ordem (precisa ser sequencial)
    for( pivo_col = 0; pivo_col < order; pivo_col++){
    	
		local_pivo.val = -1;
    	local_pivo.ind = -1;

		// Achar pivo local entre linhas nao avaliadas
    	for( j = 0; j < num_rows; j++){
    		if(row_status[j] == 0){
    			if( fabs(recvbuf[j*order + pivo_col]) > local_pivo.val ){
    				local_pivo.val = fabs(recvbuf[j*order + pivo_col]);
    				local_pivo.ind = ind_first_row + j;
    			}
    		}
    	}

    	// Achar pivo global -> reduce (MAXLOC)
    	MPI_Reduce( &local_pivo, &pivo_reduce, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD ); 

    	if(my_rank == 0){ // Armazena posicao para imprimir na ordem correta
    		pivo_order[pivo_col] = pivo_reduce.ind;	
    	}
    	

    	// Broadcast do pivo global
    	MPI_Bcast(&pivo_reduce, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);

    	// Atualizar linha do pivo
    	int local_pivo_row = pivo_reduce.ind - ind_first_row;
    	int aux = pivo_reduce.ind / ((order/num_proc) + 1);
		int pivo_rank = aux < remainder? aux : remainder + ( (pivo_reduce.ind - ((order/num_proc) + 1) * remainder) / (order/num_proc));

		// Apenas no processo que contem o pivo global
    	if(pivo_reduce.ind >= ind_first_row && pivo_reduce.ind <= ind_first_row+num_rows){
			// Atualiza o status da linha
    		row_status[local_pivo_row] = 1;

			// Copia a linha do pivo para o buffer
			#pragma omp parallel num_threads(NUM_THREADS) shared (order, local_pivo_row, pivo_row, recvbuf)
			{
				#pragma omp for
	    		for( j = 0; j < order; j++){
		    		pivo_row[j] = recvbuf[local_pivo_row * order + j];
		    	}
			}
	    	pivo_row_v = recvbuf_v[local_pivo_row];
			// Atualiza a linha no buffer
	    	#pragma omp parallel num_threads(NUM_THREADS) shared (order, pivo_row, pivo_col)
			{
	    		#pragma omp for
	    			for (int k = pivo_col+1; k < order; k++){
	    				pivo_row[k] = pivo_row[k]/pivo_row[pivo_col];
	    			}

				pivo_row_v /= pivo_row[pivo_col];
   	    		pivo_row[pivo_col] = 1;
			}
    	}


    	// Broadcast da linha atualizada
    	MPI_Bcast(pivo_row, order, MPI_DOUBLE, pivo_rank, MPI_COMM_WORLD);
    	MPI_Bcast(&pivo_row_v, 1, MPI_DOUBLE, pivo_rank, MPI_COMM_WORLD);

    	// Atualizar demais linhas
		for(j = 0; j < num_rows; j++){

			// Linhas diferentes do pivo fazem o calculo do novo valor
    		if( ind_first_row + j != pivo_reduce.ind ){	
    			double factor = recvbuf[j * order + pivo_col];

				recvbuf_v[j] = recvbuf_v[j] - (factor * pivo_row_v);
    			#pragma omp parallel num_threads(NUM_THREADS) shared (num_rows, recvbuf, pivo_row, pivo_col, factor, j, order)
				{
		    		#pragma omp for
		    		for (k = pivo_col+1; k < order; k++){
		    			recvbuf[j * order + k] = recvbuf[j * order + k] - (factor * pivo_row[k]);
		    		}
		    		
		    		recvbuf[j * order + pivo_col] = 0;
		    	}
			// A linha do pivo apenas copia o conteudo da linha atualizada
	    	} else {
    			#pragma omp parallel num_threads(NUM_THREADS) shared (recvbuf, pivo_row, pivo_col, j, order)
				{
					#pragma omp for
		    		for(k = pivo_col; k < order; k++){
		    			recvbuf[j * order + k] = pivo_row[k];
		    		}
				}

	    		recvbuf_v[local_pivo_row] = pivo_row_v;
	    	}
    	}
    }

	

    // Coleta (Gather) as informacoes do vetor resultado
    MPI_Gatherv(recvbuf_v, num_rows, MPI_DOUBLE, vector, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Escreve o resultado
    if(my_rank == 0){
		/*
		// Escreve o tempo de execucao calculado
		time_f = omp_get_wtime();
	    time_f -= time_i;
		printf("tempo: %lf\n", time_f);
		*/

		FILE *fp = fopen("resultado.txt", "w");
		for(j = 0; j < order; j++){
            fprintf(fp, "%.3lf\n", vector[pivo_order[j]]);
        }
		fclose(fp);
    }

	// Liberacao de memoria alocada
    free(recvbuf);
	free(row_status);
	free(pivo_row);
	free(recvbuf_v);
	if(my_rank == 0){
	    free(vector);
	    free(matrix);
	    free(sendcounts);
    	free(displs);
    	free(pivo_order);
	}

    MPI_Finalize();
    return 0;
}

char * ReadLine(FILE * fp, int * dim){
    char * line = NULL;
    char c;
    int count =0;
    
    do {    
        c = fgetc(fp);
        
        if(c != '\n' && c != EOF){
            line = (char*)realloc(line, sizeof(char)*(count+1));
            line[count++] = c;
        }else if(count > 0){
            line = (char*)realloc(line,sizeof(char)*(count+1));
            line[count] = '\0';   
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

    char *line =NULL;
    double *matrix = NULL;
    int i ,k = 0,order = 0, length=0;
    FILE *fp = fopen(filename, "r");
    char *p, *q;
    
    line = ReadLine(fp, &length);
    order = getOrder(line, length);
    *dim = order;
    matrix = (double*)calloc(sizeof(double), order * order);

    do {
       p = &line[0];
       i = 0;
       
       while(i < order){
          matrix[(k * order) + i] = strtod(p, &q);
		  i++;
          q++;
          p = q;
       }

       k++;
       free(line);
       line = ReadLine(fp, &length);
    }while(line != NULL); 

    fclose(fp);
    return matrix;
}

double * getVector (char * filename, int size){
        FILE * fp = fopen(filename, "r");
		if(fp == NULL){
			printf("invalide filename\n");
			return NULL;
		}
        int i = 0;
        double * vector = (double*)calloc(sizeof(double), size);

        while(i < size){
                fscanf(fp, "%lf", &vector[i]);
                i++;
        }    
        fclose(fp);
        return vector;
}
