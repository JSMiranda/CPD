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
	char *cstr1, *cstr2;

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
	
	MPI_Bcast(&nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	cstr1 = (char *) malloc(nRows*sizeof(char));
	strcpy(cstr1, sequence1.c_str());
	cstr2 = (char *) malloc(nCols*sizeof(char));
	strcpy(cstr2, sequence2.c_str());
	MPI_Bcast(cstr1, nRows, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(cstr2, nCols, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	chunkLenght = (int) floor(1.0*nCols/p/EFFICIENCY_FACTOR);
	chunksPerRow = (int) ceil(nCols*1.0/chunkLenght);
	chunksPerCol = (int) ceil(nRows*1.0/chunkLenght);
	chunksPerProcessor = chunksPerRow*chunksPerCol/p;
	
	if(id == 0){
		std::cout << "chunkLenght: " << chunkLenght << std::endl
		<< "chunkPerRow: " << chunksPerRow << std::endl
		<< "chunksPerCol: " << chunksPerCol << std::endl
		<< "chunksPerProcessor: " << chunksPerProcessor << std::endl;
		//return 0;
	}
	
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
	//int countChunks = 0;
	
	for (int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
	    //std::cout << ++countChunks << std::endl;
        bool isFirstLineLCSmatrix = (id == 0 && chunk < chunksPerRow); 
		if(!isFirstLineLCSmatrix) {
			std::cout << "(" << id << ")" << "Waiting to rcv from processor " << ((id == 0) ? (p-1) : (id - 1)) << " || Tag " << chunk%chunksPerRow << "\n";
            MPI_Recv(chunkReceiver[chunk], chunkLenght, MPI_INT, (id == 0) ? p-1 : id - 1, chunk%chunksPerRow, MPI_COMM_WORLD, &status);
			//std::cout << "(" << id << ")" << "Rcv!\n";
        }
		
		// para processar i = 0
		const bool firstRowLCSMat = id == 0 && chunk < chunksPerRow;
		if(firstRowLCSMat) {
			for(int j = 0 ; j < chunkLenght ; j++) {
				chunkArray[chunk][0][j] = 0;
			}
        // se for outra, utilizar o chunkReceiver
		} else {
			for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < nRows ; j++) {
				if (cstr1[0+chunk*chunkLenght-1] == cstr2[id + p*(int) floor(chunk / chunksPerRow)]) {
				    chunkArray[chunk][0][j] = chunkReceiver[chunk][j-1] + cost(j);
				} else if (chunkReceiver[chunk][j] >= chunkArray[chunk][0][j-1]) {
					chunkArray[chunk][0][j] = chunkReceiver[chunk][j];
				} else {
					chunkArray[chunk][0][j] = chunkArray[chunk][0][j-1];
				}
			}
		}

		// para processar j = 0
		const bool firstColLCSMat = chunk%chunksPerRow == 0;
        if(firstColLCSMat) {
			for(int i = 0; i < chunkLenght; i++) {
				chunkArray[chunk][i][0] = 0;
			}
        } else {
        // se for outra, utilizar o chunkAnterior
			for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < nCols ; i++) {
				if (cstr1[i+chunk*chunkLenght-1] == cstr2[id + p*(int) floor(chunk / chunksPerRow)]) {
			        chunkArray[chunk][i][0] = chunkArray[chunk-1][i-1][chunkLenght-1] + cost(i);
				} else if (chunkArray[chunk][i-1][0] >= chunkArray[chunk-1][i][chunkLenght-1]) {
					chunkArray[chunk][i][0] = chunkArray[chunk][i-1][0];
				} else {
					chunkArray[chunk][i][0] = chunkArray[chunk-1][i][chunkLenght-1];
				}
			}
		}

        // caso especial entre os especiais: i=0 && j=0
		if (firstRowLCSMat || firstColLCSMat) {
		    // do nothing
		} else if (cstr1[0+chunk*chunkLenght-1] == cstr2[id + p*(int) floor(chunk / chunksPerRow)]) {
		    chunkArray[chunk][0][0] = chunkReceiver[chunk-1][chunkLenght-1] + cost(chunk);
		} else if (chunkReceiver[chunk][0] >= chunkArray[chunk-1][0][chunkLenght-1]) {
			chunkArray[chunk][0][0] = chunkReceiver[chunk][0];
		} else {
			chunkArray[chunk][0][0] = chunkArray[chunk-1][0][chunkLenght-1];
		}

		if (chunk == 41)
			std::cout << id << "HI!\n";
		// Restantes casos	
		for(int i = 1; i < chunkLenght && i+chunk*chunkLenght-1 < nCols ; i++) {
		    //std::cout << chunk << ":" << i << "   nCols:" << nCols << std::endl ;
			for(int j = 1; j < chunkLenght && id + p*floor(chunk / chunksPerRow) < nRows ; j++) {
				//std::cout << "JOTA " << chunk << ":" << j << std::endl;
				if (cstr1[i+chunk*chunkLenght-1] == cstr2[id + p*(int) floor(chunk / chunksPerRow)]) {
		            chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j-1] + cost(i);
				} else if (chunkArray[chunk][i-1][j] >= chunkArray[chunk][i][j-1]) {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j];
				} else {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i][j-1];
				}
			}
		}


		const bool isLastProcess = (id == chunksPerCol % p);
		const bool isLastRowOfChunks = (chunk >= chunksPerProcessor-chunksPerRow);
        const bool isLastRowLCSmatrix = isLastProcess && isLastRowOfChunks;
        if (!isLastRowLCSmatrix) {
			std::cout << "(" << id << ")" << "Sending to processor " << (id+1)%p << " || Tag " << chunk%chunksPerRow << "\n";
            MPI_Isend(chunkArray[chunk], chunkLenght, MPI_INT, (id+1)%p, chunk%chunksPerRow, MPI_COMM_WORLD, &req);
        } else {
			std::cout << "(" << id << ")" << "HI!\n";
		}		
    }
	
	std::cout << "Algum processo parou!\n";

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
			printf("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
		}
	}
    MPI_Finalize();
    return 0;
}
