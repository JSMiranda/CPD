#define EFFICIENCY_FACTOR 20

#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stack>
#include <math.h>

std::string sequence1, sequence2;
int sequenceSize1, sequenceSize2;
int ***chunkArray;

short cost(int x){
   int i, n_iter = 20;
   double dcost = 0;
   for(i = 0; i < n_iter; i++)
      dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
   return (short) (dcost / n_iter + 0.1);
}

void readFile(int argc, char *argv[]) {
        /** Initializations **/
        std::string filename;
        std::ifstream file;
 
 
        /** Read file name from input **/
 
        if(argc == 2)
                filename = argv[1];
        else
                return;
 
 
        /** Read file and store sequences. **/
 
        file.open(filename.c_str());
 
        file >> sequenceSize1 >> sequenceSize2;
 
        std::getline(file, sequence1);   // ignore one line - to change
        std::getline(file, sequence1);
        std::getline(file, sequence2);
 
        file.close();    
}

int main (int argc, char *argv[]) {

    MPI_Status status;
	MPI_Request req;
    int id, p, nRows, nCols, chunkLenght, chunksPerRow, chunksPerCol, chunksPerProcessor;
    double secs;
	char* cstr1, *cstr2;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    MPI_Barrier (MPI_COMM_WORLD);
    secs = - MPI_Wtime();

	if(id == 0){
		readFile(argc, argv);
		nRows = sequenceSize1+1;
		nCols = sequenceSize2+1;
	}
	
	MPI_Bcast(&nCols, 1, MPI_INT, 1, MPI_COMM_WORLD);
    MPI_Bcast(&nRows, 1, MPI_INT, 1, MPI_COMM_WORLD);
	cstr1 = (char *) malloc(nRows*sizeof(char));
	strcpy(cstr1, sequence1.c_str());
	cstr2 = (char *) malloc(nCols*sizeof(char));
	strcpy(cstr2, sequence2.c_str());
	MPI_Bcast((char *) sequence1.c_str(), sequenceSize1+1, MPI_CHAR, 1, MPI_COMM_WORLD);
	MPI_Bcast((char *) sequence2.c_str(), sequenceSize2+1, MPI_CHAR, 1, MPI_COMM_WORLD);
	
	chunkLenght = floor(1.0*nCols/p/EFFICIENCY_FACTOR);
	chunksPerRow = ceil(nCols*1.0/chunkLenght);
	chunksPerCol = ceil(nRows*1.0/chunkLenght);
	chunksPerProcessor = chunksPerCol/p;
	
	int* array = (int *) malloc(chunksPerProcessor*chunkLenght*chunkLenght*sizeof(int));
	int** M = (int **) malloc(chunksPerProcessor*chunkLenght*sizeof(int*)); 
	chunkArray = (int ***) malloc(chunksPerProcessor*sizeof(int**));
	for(int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
		chunkArray[chunk] = &M[chunk*chunkLenght];
		for(int i = 0 ; i < chunkLenght ; i++) {
		    chunkArray[chunk][i] = &array[chunk*chunkLenght*chunkLenght + i*chunkLenght];
		}
	}
	
	int* receiverArray = (int *) malloc(chunksPerProcessor*chunkLenght*sizeof(int));
	int** chunkReceiver = (int **) malloc(chunksPerProcessor*sizeof(int*));
	for(int chunk = 0; chunk < chunksPerProcessor; chunk++) {
		chunkReceiver[chunk] = &receiverArray[chunk*chunkLenght];
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
	
	for (int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
        bool isFirstLineLCSmatrix = (id == 0 && chunk < chunksPerRow); 
		if(!isFirstLineLCSmatrix) {
            MPI_Recv(chunkReceiver[chunk], chunkLenght, MPI_INT, id, id*chunk, MPI_COMM_WORLD, &status);
        }
		
		// para processar i = 0
        // se for a primeira linha da matriz LCS final, são zeros
		if(id == 0 && chunk < chunksPerRow) {
			for(int j = 0 ; j < chunkLenght ; j++) {
				chunkArray[chunk][0][j] = 0;
			}
        // se for outra, utilizar o chunkReceiver
		} else {
			for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < nRows ; j++) {
				if (cstr1[0+chunk*chunkLenght-1] == cstr2[id + p*floor(chunk / chunksPerRow)]) {
				    chunkArray[chunk][0][j] = chunkReceiver[chunk][j-1] + cost(j);
				} else if (chunkReceiver[chunk][j] >= chunkArray[chunk][0][j-1]) {
					chunkArray[chunk][0][j] = chunkReceiver[chunk][j];
				} else {
					chunkArray[chunk][0][j] = chunkArray[chunk][0][j-1];
				}
			}
		}

		// para processar j = 0
        // se for a primeira coluna da matriz LCS final, são zeros
        if(chunk%chunksPerRow == 0) {
			for(int i = 0; i < chunkLenght; i++) {
				chunkArray[chunk][i][0] = 0;
			}
        } else {
        // se for outra, utilizar o chunkAnterior
			for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < nCols ; i++) {
				if (cstr1[i+chunk*chunkLenght-1] == cstr2[id + p*floor(chunk / chunksPerRow)]) {
			        chunkArray[chunk][i][0] = chunkArray[chunk-1][i-1][chunkLenght-1] + cost(i);
				} else if (chunkArray[chunk][i-1][0] >= chunkArray[chunk-1][i][chunkLenght-1]) {
					chunkArray[chunk][i][0] = chunkArray[chunk][i-1][0];
				} else {
					chunkArray[chunk][i][0] = chunkArray[chunk-1][i][chunkLenght-1];
				}
			}
		}

        // caso especial entre os especiais: i=0 && j=0
		if (cstr1[0+chunk*chunkLenght-1] == cstr2[id + p*floor(chunk / chunksPerRow)]) {
		    chunkArray[chunk][0][0] = chunkReceiver[chunk-1][chunkLenght-1] + cost(chunk);
		} else if (chunkReceiver[chunk][0] >= chunkArray[chunk-1][0][chunkLenght-1]) {
			chunkArray[chunk][0][0] = chunkReceiver[chunk][0];
		} else {
			chunkArray[chunk][0][0] = chunkArray[chunk-1][0][chunkLenght-1];
		}

		// Restantes casos	
		for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < nCols ; i++) {
			for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < nRows ; j++) {
				if (cstr1[i+chunk*chunkLenght-1] == cstr2[id + p*floor(chunk / chunksPerRow)]) {
		            chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j-1] + cost(i);
				} else if (chunkArray[chunk][i-1][j] >= chunkArray[chunk][i][j-1]) {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j];
				} else {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i][j-1];
				}
			}
		}


        bool isLastLineLCSmatrix = chunksPerCol % p;
        if (!isLastLineLCSmatrix) {
            MPI_Isend(chunkArray[chunk], chunkLenght, MPI_INT, (id+1)%p, ((id+1)%p)*chunk, MPI_COMM_WORLD, &req);
        }	
    }

    MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();

    if(id == 0){
       printf("Processes = %d, Time = %12.6f sec,\n", p, secs);
	   printf("%d\n", chunksPerProcessor);
		printf("(p0) String1: %s\n", cstr1);
	   for (int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
		   for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < nCols ; i++) {
				for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < nRows ; j++) {
					printf("%d ", chunkArray[chunk][i][j]);
				}
				printf("\n");
			}
		}
    } else {
		printf("String1: %s\n", cstr1);
	}
    MPI_Finalize();
    return 0;
}
