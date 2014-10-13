#include <math.h>
#include <stdio.h>
#include <omp.h>

#define ENDOFLOOP 2000000
#define TRUE 1
#define FALSE 0

//int x;
//#pragma omp threadprivate(x)

short cost(int x)
{
   int i, n_iter = 20;
   double dcost = 0;
   for(i = 0; i < n_iter; i++)
      dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
   return (short) (dcost / n_iter + 0.1);
}

/*
 * To run with OpenMP, run ./test p
 * Otherwise will run sequentially
 */
void main(int argc, char *argv[]) {
   int i;
   if(argc == 2 && !strcmp(argv[1], "p")) {
      #pragma omp parallel
      {
          //x = omp_get_thread_num();
          #pragma omp for
          for(i = 0; i < ENDOFLOOP ; i++) {
	         cost(i);
	         //printf("Running thread %d\n", tid);
          }
       }
   } else {
	   for(i = 0; i < ENDOFLOOP ; i++) {
	      cost(i);
	   }
   }
}
