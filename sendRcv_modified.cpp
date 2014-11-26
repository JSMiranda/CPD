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
    int id, p, nRows, nCols, chunkLenght, chunksPerRow, chunksPerCol, chunksPerProcessor;
    double secs;

    // FIXME: verificar isto
    MPI_Init (&argc, &argv);

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
	MPI_Bcast(sequence1.c_str(), sequenceSize1+1, MPI_CHAR, 1, MPI_COMM_WORLD);
	MPI_Bcast(sequence2.c_str(), sequenceSize2+1, MPI_CHAR, 1, MPI_COMM_WORLD);
	
	chunkLenght = floor(1.0*nCols/p/EFFICIENCY_FACTOR);
	chunksPerRow = ceil(nCols*1.0/chunkLenght);
	chunksPerCol = ceil(nRows*1.0/chunkLenght);
	chunksPerProcessor = chunksPerCol/p;
	
	int* array = (int *) malloc(chunksPerProcessor*chunkLenght*chunkLenght*sizeof(int));
	int** M = (int **) malloc(chunksPerProcessor*chunkLenght*sizeof(int*)); 
	chunkArray = (int ***) malloc(chunksPerProcessor*sizeof(int**));
	for(int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
		chunkArray[chunk] = &M[chunk*chunkLength];
		for(int i = 0 ; i < chunksLenght ; i++) {
		    chunkArray[chunk][i] = &array[chunk*chunkLenght*chunkLenght + i*chunkLenght];
		}
	}
	
	int* receiver_array = (int *) malloc(chunksPerProcessor*chunkLenght*sizeof(int));
	int** chunkReceiver = (int **) malloc(chunksPerProcessor*sizeof(int*));
	for(int chunk = 0; chunk < chunksPerProcessor; chunk++) {
		chunkReceiver[chunk] = &receiver_array[chunk*chunkLenght];
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
	
	for (int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
        bool isFirstLineLCSmatrix = (id == 0 && chunk < chunksPerRow); 
		if(!isFirstLineLCSmatrix) {
            // TODO: recebe a ultima linha do chunk do processo anterior
        }
		
		// para processar i = 0
        // se for a primeira linha da matriz LCS final, são zeros
		if(id == 0 && chunk < chunksPerRow) {
			for(int j = 0 ; j < chunkLenght) {
				chunkArray[chunk][0][j] = 0;
			}
        // se for outra, utilizar o receiverArray
		} else {
			for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < rows ; j++) {
				if (sequence1.at(0+chunk*chunkLenght-1) == sequence2.at(id + p*floor(chunk / chunksPerRow))) {
				    chunkArray[chunk][0][j] = receiverArray[chunk][j-1] + cost(j);
				} else if (receiverArray[chunk][j] >= chunkArray[chunk][0][j-1]) {
					chunkArray[chunk][0][j] = receiverArray[chunk][j];
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
			for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < cols ; i++) {
				if (sequence1.at(i+chunk*chunkLenght-1) == sequence2.at(id + p*floor(chunk / chunksPerRow))) {
			        chunkArray[chunkID][i][0] = chunkArray[chunk-1][i-1][chunkSize-1] + cost(i);
				} else if (chunkArray[chunkID][i-1][0] >= chunkArray[chunk-1][i][chunkSize-1]) {
					chunkArray[chunkID][i][0] = chunkArray[chunk][i-1][0];
				} else {
					chunkArray[chunkID][i][0] = chunkArray[chunk-1][i][chunkSize-1];
				}
			}
		}

        // caso especial entre os especiais: i=0 && j=0
		if (sequence1.at(0+chunk*chunkLenght-1) == sequence2.at(id + p*floor(chunk / chunksPerRow))) {
		    chunkArray[chunk][0][0] = receiverArray[chunk-1][chunkSize-1] + cost(chunk);
		} else if (receiverArray[chunk][j] >= chunkArray[chunk-1][0][chunkSize-1]) {
			chunkArray[chunk][0][0] = receiverArray[chunk][j];
		} else {
			chunkArray[chunk][0][0] = chunkArray[chunk-1][0][chunkSize-1];
		}

		// Restantes casos	
		for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < cols ; i++) {
			for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < rows ; j++) {
				if (sequence1.at(i+chunk*chunkLenght-1) == sequence2.at(id + p*floor(chunk / chunksPerRow))) {
		            chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j-1] + cost(i);
				} else if (chunkArray[chunk][i-1][j] >= chunkArray[chunk][i][j-1]) {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j];
				} else {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i][j-1];
				}
			}
		}


        bool isLastLineLCSmatrix = (ultimoProcesso = chunksPerCol % p);
        if (!isLastLineLCSmatrix) {
		    // envia ultima linha do chunk ao processo seguinte
        }	
    }

    MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();

    if(id == 0){
       printf("Processes = %d, Time = %12.6f sec,\n", p, secs);
    }
    MPI_Finalize();
    return 0;
}
