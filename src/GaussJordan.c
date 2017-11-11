/*
// para compilar: mpicc GaussJordan.c -o GaussJordan -Wall -lm -fopenmp
// para rodar: mpirun -np 2 GaussJordan
*/
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <omp.h>

int getOrder(char * line, int length);
char * ReadLine(FILE *input, int * dim);
double ** getMatrix(char * filename, int * dim );
double * getVector (char * filename, int size);
void destroyMatrix(double ** vector, int order);

int main(int argc, char *argv[]){
    int *result;
    int order;
    double ** matrix = NULL;
    double * vector=NULL;
    int my_rank, num_proc;      
    double * recvbuf;
    int * sendcounts = NULL, * displs = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

/*A-ler arquivos: matrix e vetor*/
    
    if(my_rank == 0){
        printf("sou rank 0\n");
        matrix = getMatrix("matriz.txt", &order);
        vector = getVector("vetor.txt", order);
    }
    	
    
/*B-scatter da matriz/vetor (varias linhas por processo -> COMO FAZER A DIVISÃO?)*/

    
    
    int num_elements = (order / num_proc) * order;
    int rest = (order % num_proc);
    
    recvbuf = (double *)malloc(sizeof(double) * (num_elements + order));
    
    if(my_rank == 0){
        sendcounts = (int *)malloc(sizeof(int) * num_proc);
        displs = (int *)malloc(sizeof(int) * num_proc);
            
        int sum = 0;
        
        for(int i = 0; i < num_proc; i++){
            sendcounts[i] = num_elements;
            
            if(rest > 0){
                sendcounts[i] += order;
                rest--;
            }
            
            displs[i] = sum;
            sum += sendcounts[i];
        }
    }
    
    MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, &recvbuf, num_elements + order, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int aux = order - my_rank;
    if(aux < 0) aux = 0;
    
    for(int i = 0; i < num_elements + aux; i++)
        printf("meu rank: %d", my_rank);

/*C-Gauss jordan*/

    /*Passo 1: encontrar pivot */
    
    /*Passo 2: trocar linhas*/

    /*Passo 3: dividir elementos da linha k pelo pivot - n tarefas */

    /*Passo 4: somar linha k com os valores das demais linhas */
    
    
    /*free(recvbuf);
    if(sendcounts != NULL)
    	free(sendcounts);
    if(displs != NULL)
    	free(displs);
    */
    destroyMatrix(matrix, order);
    free(vector);

    MPI_Finalize();
    return 0;
}


void destroyMatrix(double ** matrix, int order){
    int i;
    for (i=0; i<order; i++){
        
        free(matrix[i]);
    
    }
    free(matrix);

}

char * ReadLine(FILE * fp, int * dim){
    char * line = NULL;
    char c;
    int count =0;
    
    do {    
        c = fgetc(fp);
        
        if(c!='\n' && c!=EOF){
            line = (char*)realloc(line, sizeof(char)*( count+1));
            line [count++] = c;
             
        }else if(count>0){
            line = (char*)realloc(line,sizeof(char)*(count+1));
            line [count++] = '\0';   
        
        }
    
    }while(c!=EOF && c!='\n');
 
    *dim = count;
    return line;
}


int getOrder(char * line, int length){
    int count =1;
    int i=0;

    while(i<length){
        if(line[i]==' ') count ++;
        i++;
    }
    

    return count;
}

double ** getMatrix(char * filename, int * dim){

    char * line =NULL;

    double ** matrix = NULL;

    int i ,k = 0,order = 0, length=0;

    FILE * fp = fopen(filename, "r");

    char * p, *q;

    
    
    line = ReadLine(fp, &length);
    
    order = getOrder (line, length);

    *dim = order;

    matrix = (double**)calloc(sizeof(double*), order);
      
    
        
    do {
               
       matrix[k] = (double*)calloc(sizeof(double), order);
         
       p=&line[0]; 
      
       i=0;
       
       while(i<order){
    
            matrix[k][i++]=strtod(p, &q);
            q++;
        p = q;
                    
       }
       
       k++; 
       
       free(line);   
      
       line = ReadLine(fp, &order);
       
    
    }while(line!=NULL); 


    fclose(fp);
    return matrix;
}
double * getVector (char * filename, int size){
        FILE * fp = fopen(filename, "r");
        int i=0;
        double * vector = (double*)calloc(sizeof(double), size);

        while(i<size){
                fscanf(fp, "%lf", &vector[i]);
                i++;
        }    
        fclose(fp);
        return vector;
}

