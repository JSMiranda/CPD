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
		
		chunkLenght = floor(1.0*nCols/p/EFFICIENCY_FACTOR);
        chunksPerRow = ceil(nCols*1.0/chunkLenght);
        chunksPerCol = ceil(nRows*1.0/chunkLenght);
		chunksPerProcessor = chunksPerCol/p;
	}
	
	MPI_Bcast(&nCols, 1, MPI_INT, 1, MPI_COMM_WORLD);
    MPI_Bcast(&chunkLenght, 1, MPI_INT, 1, MPI_COMM_WORLD);
    MPI_Bcast(&chunksPerProcessor, 1, MPI_INT, 1, MPI_COMM_WORLD);

	MPI_Bcast(sequence1.c_str(), sequenceSize1+1, MPI_CHAR, 1, MPI_COMM_WORLD);
	MPI_Bcast(sequence2.c_str(), sequenceSize2+1, MPI_CHAR, 1, MPI_COMM_WORLD);
	
	int* array = (int *) malloc(chunksPerProcessor*chunkLenght*chunkLenght*sizeof(int));
	int*** chunkArray = (int ***) malloc(chunksPerProcessor*sizeof(int**));
	for(int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
		chunkArray[chunk] = (int **) malloc(chunkLenght*sizeof(int*));
		for(int i = 0 ; i < chunksLenght ; i++) {
		    chunkArray[chunk][i] = &array[chunk*chunkLenght*chunkLenght + i*chunkLenght];
		}
	}
	
	int* receiver_array = (int *) malloc(chunksPerProcessor*chunkLenght*sizeof(int));
	int** chunkReceiver = (int **) malloc(chunksPerProcessor*sizeof(int*));
	for(int chunk = 0; chunk < chunksPerProcessor; chunk++) {
		chunkReceiver[chunk] = &receiver_array[chunk*chunkLenght];
	}
	
	// first row = 0
	if(id == 0) {
		for(chunk = 0; chunk < chunksPerRow; chunk++) {
			for(i = 0; i < chunkLenght; i++) {
				chunkArray[chunk][0][i] = 0;
			}
		}
	}

    // first column = 0
	for(chunk = 0; chunk < chunksPerProcessor; chunk += chunksPerRow) {
		for(i = 0; i < chunkLenght; i++) {
			chunkArray[chunk][i][0] = 0;
		}
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
	
	// for each chunk
		// exceto na primeira linha, recebe a ultima linha do chunk do processo anterior
		// Processa chunk
		// exceto na ultima linha envia ultima linha do chunk ao processo seguinte
		
		
	/*
	 * ultimoProcesso = chunksPerCol % p;
	 */
	
	
//	else{
//	    MPI_Recv(&i, 1, MPI_INT, id-1, i, MPI_COMM_WORLD, &status);
//	    MPI_Send(&i, 1, MPI_INT, (id+1)%p, i, MPI_COMM_WORLD);
//	}


    MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();

    if(id == 0){
       printf("Processes = %d, Time = %12.6f sec,\n", p, secs);
    }
    MPI_Finalize();
    return 0;
}
