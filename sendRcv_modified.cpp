#define EFFICIENCY_FACTOR 10

#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stack>
#include <math.h>
#include <unistd.h>
#include <algorithm>

std::string sequence1, sequence2;
int sequenceSize1, sequenceSize2;
int ***chunkArray;
int** chunkReceiver;
MPI_Status status;
MPI_Request req;
int id, p, nRows, nCols, chunkLenght, chunksPerRow, chunksPerCol, chunksPerProcessor, chunkRowsPerProcessor;
double secs;
char *cstr1, *cstr2;

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

int getLCSrow(int chunk, int offset) {
	if(chunk/chunksPerRow > chunksPerCol/p) {
		return p*(chunkRowsPerProcessor-1)*chunkLenght + id*chunkLenght + offset;
    } else {
		return id*chunkLenght + p*(chunk/chunksPerRow)*chunkLenght + offset;
    }
}

void processChunks() {
	for (int chunk = 0 ; chunk < chunksPerProcessor ; chunk++) {
		const bool firstRowLCSMat = (id == 0 && chunk < chunksPerRow);
		const bool firstColLCSMat = (chunk%chunksPerRow == 0);
		
		if(!firstRowLCSMat) {
			//std::cout << "(" << id << ")" << "Waiting to rcv from processor " << ((id == 0) ? (p-1) : (id - 1)) << " || Tag " << chunk%chunksPerRow << "\n";
            MPI_Recv(chunkReceiver[chunk], chunkLenght, MPI_INT, (id == 0) ? p-1 : id - 1, chunk%chunksPerRow, MPI_COMM_WORLD, &status);
			//std::cout << "(" << id << ")" << "Rcv! CHUNK " << chunk << "   ";
			//for(int i = 0 ; i < chunkLenght ; i++) {
			//    std::cout << chunkReceiver[chunk][i] << " ";
			//}
			//std::cout << std::endl;
        }
		
		// caso especial entre os especiais: i=0 && j=0
		if (firstRowLCSMat || firstColLCSMat) {
		    // do nothing
		} else if (cstr1[getLCSrow(chunk, 0) - 1] == cstr2[(chunk%chunksPerRow)*chunkLenght - 1]) {
		    chunkArray[chunk][0][0] = chunkReceiver[chunk-1][chunkLenght-1] + cost(chunk);
		} else if (chunkReceiver[chunk][0] >= chunkArray[chunk-1][0][chunkLenght-1]) {
			chunkArray[chunk][0][0] = chunkReceiver[chunk][0];
		} else {
			chunkArray[chunk][0][0] = chunkArray[chunk-1][0][chunkLenght-1];
		}
		
		// para processar i = 0
		if(firstRowLCSMat) {
			for(int j = 0 ; j < chunkLenght ; j++) {
				chunkArray[chunk][0][j] = 0;
			}
		} else {
			for(int j = 1; j < chunkLenght && j+(chunk%chunksPerRow)*chunkLenght < nCols ; j++) {
			
				/////////////////////
				if ( j+(chunk%chunksPerRow)*chunkLenght >= nCols) {
					}
				/////////////////////
			
			
				if (cstr1[getLCSrow(chunk, 0) - 1] == cstr2[j+(chunk%chunksPerRow)*chunkLenght - 1]) {
				    chunkArray[chunk][0][j] = chunkReceiver[chunk][j-1] + cost(j);
				} else if (chunkReceiver[chunk][j] >= chunkArray[chunk][0][j-1]) {
					chunkArray[chunk][0][j] = chunkReceiver[chunk][j];
				} else {
					chunkArray[chunk][0][j] = chunkArray[chunk][0][j-1];
				}
			}
		}

		// para processar j = 0
        if(firstColLCSMat) {
			for(int i = 0; i < chunkLenght; i++) {
				chunkArray[chunk][i][0] = 0;
			}
        } else {
			for(int i = 1; i < chunkLenght && getLCSrow(chunk, i) < nRows ; i++) {
			
			/////////////////////
				if ( getLCSrow(chunk, i) >= nRows ) {
					}
				/////////////////////
			
			
				if (cstr1[getLCSrow(chunk, i) - 1] == cstr2[(chunk%chunksPerRow)*chunkLenght - 1]) {
			        chunkArray[chunk][i][0] = chunkArray[chunk-1][i-1][chunkLenght-1] + cost(i);
				} else if (chunkArray[chunk][i-1][0] >= chunkArray[chunk-1][i][chunkLenght-1]) {
					chunkArray[chunk][i][0] = chunkArray[chunk][i-1][0];
				} else {
					chunkArray[chunk][i][0] = chunkArray[chunk-1][i][chunkLenght-1];
				}
			}
		}

		// Restantes casos
		for(int i = 1; i < chunkLenght && getLCSrow(chunk, i) < nRows ; i++) {
			for(int j = 1; j < chunkLenght && j+(chunk%chunksPerRow)*chunkLenght < nCols ; j++) {
				if (cstr1[getLCSrow(chunk, i) - 1] == cstr2[j+(chunk%chunksPerRow)*chunkLenght - 1]) {
		            chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j-1] + cost(i);
					/*
						std::cout << "Found match at " << getLCSrow(chunk, i) << ", " << j+(chunk%chunksPerRow)*chunkLenght << " : "
						<< cstr1[getLCSrow(chunk, i) - 1] << std::endl << "Chunk " << chunk << " | i " << i << " | j " << j << " Value: " << chunkArray[chunk][i][j] << std::endl;
					
					*/
				} else if (chunkArray[chunk][i-1][j] >= chunkArray[chunk][i][j-1]) {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i-1][j];
				} else {
				    chunkArray[chunk][i][j] = chunkArray[chunk][i][j-1];
				}
			}
		}


		const bool isLastProcess = (id == (chunksPerCol%p - 1 + p)%p);
		const bool isLastRowOfChunks = (chunk >= chunksPerProcessor-chunksPerRow);
        const bool isLastRowLCSmatrix = isLastProcess && isLastRowOfChunks;
        if (!isLastRowLCSmatrix) {
			//std::cout << "(" << id << ")" << "Sending to processor " << (id+1)%p << " || Tag " << chunk%chunksPerRow << "\n";
            MPI_Isend(chunkArray[chunk][chunkLenght-1], chunkLenght, MPI_INT, (id+1)%p, chunk%chunksPerRow, MPI_COMM_WORLD, &req);
			//std::cout << "(" << id << ")" << "Send! CHUNK " << chunk << "   ";
			//if(chunkArray[chunk][chunkLenght-1][0] == -1) {
			//	for(int i = 0 ; i < chunkLenght ; i++) {
			//		for(int j = 0 ; j < chunkLenght ; j++) {
			//			std::cout << chunkArray[chunk][i][j] << " ";
			//		}
			//	}
			//}
			//std::cout << std::endl;
		}
    }
}

void printLCS() {

	const int maxSize = std::min(nRows, nCols);
	char* result = (char *) malloc(maxSize*sizeof(char));
	
	int i, j, chunk=chunksPerProcessor-1;
	int currentSize = 0;
	
	MPI_Request req1, req2, req3, req4, req5;
	const bool isLastProcess = (id == (chunksPerCol%p - 1 + p)%p);
	bool chunkRowBeeingProcessed = false;
	
	if(isLastProcess) {
		printf("%d\n", chunkArray[chunk][nRows - getLCSrow(chunksPerProcessor-1, 0) - 1][nCols - 1 - (chunksPerRow-1)*chunkLenght]);
		i = nRows - getLCSrow(chunksPerProcessor-1, 0) - 1;
		j = nCols - 1 - (chunksPerRow-1)*chunkLenght;
		chunkRowBeeingProcessed = true;
	}
	
	bool jobDone = false;
	do {
		const bool isLastRowOfChunks = (chunk >= chunksPerProcessor-chunksPerRow);
		const bool isLastRowLCSmatrix = isLastProcess && isLastRowOfChunks;
		if(!chunkRowBeeingProcessed) {
			printf("Waiting for receive. PID:%d\n", id);
			MPI_Irecv(&i, 1, MPI_INT, (id+1)%p, 0, MPI_COMM_WORLD, &req1);
			MPI_Irecv(&j, 1, MPI_INT, (id+1)%p, 1, MPI_COMM_WORLD, &req2);
			MPI_Irecv(&currentSize, 1, MPI_INT, (id+1)%p, 3, MPI_COMM_WORLD, &req4);
			MPI_Wait(&req4, &status);
			MPI_Irecv(result, currentSize, MPI_CHAR, (id+1)%p, 2, MPI_COMM_WORLD, &req3);
			MPI_Irecv(&chunk, 1, MPI_INT, (id+1)%p, 4, MPI_COMM_WORLD, &req5);
			
			
			MPI_Wait(&req1, &status);
			MPI_Wait(&req2, &status);
			MPI_Wait(&req3, &status);
			MPI_Wait(&req5, &status);
			printf("Received. PID:%d\n", id);

		}
		
		chunkRowBeeingProcessed = true;
	
		while( i >= 1 && j >= 1 ){
			if (cstr1[getLCSrow(chunk,i)-1] == cstr2[j+(chunk%chunksPerRow)*chunkLenght-1]){
				result[currentSize] = cstr1[getLCSrow(chunk,i)-1];
				i--;
				j--;
				currentSize++;
			} else if (chunkArray[chunk][i-1][j] <= chunkArray[chunk][i][j-1]) {
				j--;
			} else {
				i--;
			}
		}
		
		const bool firstRowLCSMat = (id == 0 && chunk < chunksPerRow);
		const bool firstColLCSMat = (chunk%chunksPerRow == 0 && j==0);
		if(id == 0 && (firstRowLCSMat || firstColLCSMat)) {
			printf("CurrentSize: %d\n", currentSize);
			while(currentSize >= 0) {
				printf("%c", result[--currentSize]);
			}
			printf("\n");
			return;
		} else if(firstRowLCSMat || firstColLCSMat) {
			MPI_Isend(&i, 1, MPI_INT, (id-1+p)%p, 0, MPI_COMM_WORLD, &req1);
			MPI_Isend(&j, 1, MPI_INT, (id-1+p)%p, 1, MPI_COMM_WORLD, &req2);
			MPI_Isend(&currentSize, 1, MPI_INT, (id-1+p)%p, 3, MPI_COMM_WORLD, &req4);
			MPI_Isend(result, currentSize, MPI_CHAR, (id-1+p)%p, 2, MPI_COMM_WORLD, &req3);
			MPI_Isend(&chunk, 1, MPI_INT,(id-1+p)%p, 4, MPI_COMM_WORLD, &req5);
			jobDone=true;
		} else if(i == 0 && j == 0) {
			if (cstr1[getLCSrow(chunk,i)-1] == cstr2[j+(chunk%chunksPerRow)*chunkLenght-1]){
				result[currentSize] = cstr1[getLCSrow(chunk,i)-1];
				i=chunkLenght-1;
				j=chunkLenght-1;
				chunk--;
				if(chunk < chunksPerRow) {
					jobDone=true;
				}
				if(id == 0)
					chunk -= chunksPerRow;
				currentSize++;
				MPI_Isend(&i, 1, MPI_INT, (id-1+p)%p, 0, MPI_COMM_WORLD, &req1);
				MPI_Isend(&j, 1, MPI_INT, (id-1+p)%p, 1, MPI_COMM_WORLD, &req2);
				MPI_Isend(&currentSize, 1, MPI_INT, (id-1+p)%p, 3, MPI_COMM_WORLD, &req4);
				MPI_Isend(result, currentSize, MPI_CHAR, (id-1+p)%p, 2, MPI_COMM_WORLD, &req3);
				MPI_Isend(&chunk, 1, MPI_INT,(id-1+p)%p, 4, MPI_COMM_WORLD, &req5);
				chunkRowBeeingProcessed = false;
			} else if (chunkReceiver[chunk][j] <= chunkArray[chunk-1][i][chunkLenght-1]) {
				j=chunkLenght-1;
				chunk-=1;
			} else {
				i=chunkLenght-1;
				if(chunk < chunksPerRow) {
					jobDone=true;
				}
				if(id == 0)
					chunk -= chunksPerRow;
				MPI_Isend(&i, 1, MPI_INT, (id-1+p)%p, 0, MPI_COMM_WORLD, &req1);
				MPI_Isend(&j, 1, MPI_INT, (id-1+p)%p, 1, MPI_COMM_WORLD, &req2);
				MPI_Isend(&currentSize, 1, MPI_INT, (id-1+p)%p, 3, MPI_COMM_WORLD, &req4);
				MPI_Isend(result, currentSize, MPI_CHAR, (id-1+p)%p, 2, MPI_COMM_WORLD, &req3);
				MPI_Isend(&chunk, 1, MPI_INT,(id-1+p)%p, 4, MPI_COMM_WORLD, &req5);
				chunkRowBeeingProcessed = false;
			}
		} else if(i == 0 && j != 0) {
			if (cstr1[getLCSrow(chunk,i)-1] == cstr2[j+(chunk%chunksPerRow)*chunkLenght-1]){
				result[currentSize] = cstr1[getLCSrow(chunk,i)-1];
				i=chunkLenght-1;
				j--;
				currentSize++;
				if(chunk < chunksPerRow) {
					jobDone=true;
				}
				if(id == 0)
					chunk -= chunksPerRow;
				MPI_Isend(&i, 1, MPI_INT, (id-1+p)%p, 0, MPI_COMM_WORLD, &req1);
				MPI_Isend(&j, 1, MPI_INT, (id-1+p)%p, 1, MPI_COMM_WORLD, &req2);
				MPI_Isend(&currentSize, 1, MPI_INT, (id-1+p)%p, 3, MPI_COMM_WORLD, &req4);
				MPI_Isend(result, currentSize, MPI_CHAR, (id-1+p)%p, 2, MPI_COMM_WORLD, &req3);
				MPI_Isend(&chunk, 1, MPI_INT,(id-1+p)%p, 4, MPI_COMM_WORLD, &req5);
				chunkRowBeeingProcessed = false;
			} else if (chunkReceiver[chunk][j] <= chunkArray[chunk][i][j-1]) {
				j--;
			} else {
				i=chunkLenght-1;
				if(chunk < chunksPerRow) {
					jobDone=true;
				}
				if(id == 0)
					chunk -= chunksPerRow;
				MPI_Isend(&i, 1, MPI_INT, (id-1+p)%p, 0, MPI_COMM_WORLD, &req1);
				MPI_Isend(&j, 1, MPI_INT, (id-1+p)%p, 1, MPI_COMM_WORLD, &req2);
				MPI_Isend(&currentSize, 1, MPI_INT, (id-1+p)%p, 3, MPI_COMM_WORLD, &req4);
				MPI_Isend(result, currentSize, MPI_CHAR, (id-1+p)%p, 2, MPI_COMM_WORLD, &req3);
				MPI_Isend(&chunk, 1, MPI_INT,(id-1+p)%p, 4, MPI_COMM_WORLD, &req5);
				chunkRowBeeingProcessed = false;
			}
		} else if(i != 0 && j == 0) {
			if (cstr1[getLCSrow(chunk,i)-1] == cstr2[j+(chunk%chunksPerRow)*chunkLenght-1]){
				result[currentSize] = cstr1[getLCSrow(chunk,i)-1];
				i--;
				j=chunkLenght-1;
				chunk-=1;
				currentSize++;
			} else if (chunkArray[chunk][i-1][j] <= chunkArray[chunk-1][i][chunkLenght-1]) {
				j=chunkLenght-1;
				chunk-=1;
			} else {
				i--;
			}
		}
	} while(!jobDone);
}

int main (int argc, char *argv[]) {

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
	if(chunkLenght < 1)
		chunkLenght = 1;
	chunksPerRow = (int) ceil(1.0*nCols/chunkLenght);
	chunksPerCol = (int) ceil(1.0*nRows/chunkLenght);
	chunkRowsPerProcessor = chunksPerCol/p;

	if(id < chunksPerCol%p)
		chunkRowsPerProcessor++;
	
	chunksPerProcessor = chunkRowsPerProcessor*chunksPerRow;
	
	//Debug
	if(id == 0){
		std::cout << "chunkLenght: " << chunkLenght << std::endl
		<< "chunkPerRow: " << chunksPerRow << std::endl
		<< "chunksPerCol: " << chunksPerCol << std::endl
		<< "chunksPerProcessor: " << chunksPerProcessor << std::endl;
	}
	
	//Allocations
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
	chunkReceiver = (int **) malloc(chunksPerProcessor*sizeof(int*));
	for(int chunk = 0; chunk < chunksPerProcessor; chunk++) {
		chunkReceiver[chunk] = &receiverArray[chunk*chunkLenght];
	}
	
	// debug
	for(int c = 0 ; c < chunksPerProcessor ; c++) {
		for(int i = 0 ; i < chunkLenght ; i++) {
			for(int j = 0 ; j < chunkLenght ; j++) {
				chunkArray[c][i][j] = -1;
			}
		}
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
	processChunks();
	
	//std::cout << "Algum processo parou!\n";

    MPI_Barrier (MPI_COMM_WORLD);
	printLCS();
    secs += MPI_Wtime();

	// prints last position of LCS MATRIX
    if(id == (chunksPerCol%p - 1 + p)%p){
			int chunk = chunksPerProcessor-1;
			int i = nRows - getLCSrow(chunk, 0) - 1;
			int h = nCols - 1 - (chunksPerRow-1)*chunkLenght;
	}
	
    MPI_Finalize();
    return 0;
}
