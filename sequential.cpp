#include <iostream>
#include <fstream>
#include <string>

std::string sequence1, sequence2;
int sequenceSize1, sequenceSize2;
int **lcs_matrix;

void readFile(int argc, char *argv[]) {
        /** Initializations **/
        std::string filename;
        std::ifstream file;
 
 
        /** Read file name form input **/
 
        if(argc == 2)
                filename = argv[1];
        else
                return;
 
 
        /** Read file and store sequences. Attention to End Of File **/
 
        file.open(filename.c_str());
 
        file >> sequenceSize1 >> sequenceSize2;
 
        std::getline(file, sequence1);   // ignore one line - to chagne
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

void lcs-length() {
    //sequenceSize1;
    //sequenceSize2;
    //sequence1;
    //sequence2;
    int i, j;
    
    b = malloc(sequenceSize2 * sizeof(int));
    for(i = 0; i < sequenceSize2 ; i++) {
        b[i] = malloc(sequenceSize1 * sizeof(int));
    }

/*    
    first line and row = 0
    for each square (i,j)
        if xi == yi
            c[i][j] = c[i-1][j-1] + 1; // FIXME: cost(x)
        else if c[i-1][j] >= c[i][j-1]
            c[i][j] = c[i-1][j];
        else
            c[i][j] = c[i][j-1];
  */  
}

int main(int argc, char *argv[])
{
 
        readFile(argc, argv); 
        printTest();
 
        return 0;
}