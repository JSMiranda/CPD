#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stack>

std::string sequence1, sequence2;
int sequenceSize1, sequenceSize2;
int **lcs_matrix;

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

void printTest() {
    std::cout << sequenceSize1 << "\n";
    std::cout << sequenceSize2 << "\n";
 
    std::cout << "First Sequence:" << sequence1 << "\n\n";
    std::cout << "Second Sequence:" << sequence2 << "\n\n";
}

/*
 * Computes the lcs matrix
 */
void lcs_length() {
    int i, j;
    
    lcs_matrix = (int**) malloc((sequenceSize1+1) * sizeof(int*));
    for(i = 0; i < sequenceSize1 + 1 ; i++) {
        lcs_matrix[i] = (int*) malloc((sequenceSize2+1) * sizeof(int));
    }

    // first row = 0
	for(i = 0; i < sequenceSize2 + 1; i++) {
		lcs_matrix[0][i] = 0;
    }

    // first column = 0
	for(i = 0; i < sequenceSize1 + 1; i++) {
		lcs_matrix[i][0] = 0;
    }

    for(i = 1 ; i < sequenceSize1 + 1 ; i++) {
		for(j = 1 ; j < sequenceSize2 + 1 ; j++) {
		    if (sequence1.at(i-1) == sequence2.at(j-1)) {
		        lcs_matrix[i][j] = lcs_matrix[i-1][j-1] + 1; // FIXME: cost(x)
		    } else if (lcs_matrix[i-1][j] >= lcs_matrix[i][j-1]) {
		        lcs_matrix[i][j] = lcs_matrix[i-1][j];
			} else {
		        lcs_matrix[i][j] = lcs_matrix[i][j-1];
			}
		}
	}
  
}

void print_lcs() {
	int i, j;
    /*for(i = 0 ; i < sequenceSize1 ; i++) {
		for(j = 0 ; j < sequenceSize2 ; j++) {
			printf("%d  ", lcs_matrix[i][j]);
		}
		printf("\n");
    }*/
	printf("Size: %d\n", lcs_matrix[sequenceSize1][sequenceSize2]);


    // Print the sequence itself

    i = sequenceSize1;
    j = sequenceSize2;
    std::stack<char> subsequence;                   // stack is chosen so it inverts the subsequence

    while( i >= 1 && j >= 1 ){                             // while not at the end of one of the sequences
        if (sequence1.at(i-1) == sequence2.at(j-1)){
            subsequence.push(sequence1.at(i-1));        // saves the matching char in a stack
            i--;    // move diagonally
            j--;
        } else if (lcs_matrix[i-1][j] >= lcs_matrix[i][j-1]) {
            i--;    // move right
        } else {
            j--;    // move left
        }
    }

    // print stack

    std::cout << "Subsquence:";

    while( !subsequence.empty() )
    {
       std::cout << ' ' << subsequence.top();
       subsequence.pop();
    }

    std::cout << "\n";

}

int main(int argc, char *argv[])
{

        readFile(argc, argv); 
        printTest();
		lcs_length();
		print_lcs();
 
        return 0;
}
