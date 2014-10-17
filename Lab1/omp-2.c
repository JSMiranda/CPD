/*
DESCRIPTION:
		Parallelizing an inner loop with dependences

		for (iter=0; iter<numiter; iter++) {
			for (i=0; i<size-1; i++) {
				V[i] = f( V[i], V[i+1] );
			}
		}

**************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

#define TOTALSIZE 1000
#define NUMITER 200

/*
* DUMMY FUNCTION
*/
#define f(x,y)	((x+y)/2.0)


/* MAIN: PROCESS PARAMETERS */
int main(int argc, char *argv[]) {

  /* VARIABLES */
  int i, iter;

  /* DECLARE VECTOR AND AUX DATA STRUCTURES */
  double *V = (double *) malloc(TOTALSIZE * sizeof(double));

  /* 1. INITIALIZE VECTOR */
  for(i = 0; i < TOTALSIZE; i++) {
    V[i]= 0.0 + i;
  }

  /* 2. ITERATIONS LOOP */
  for(iter = 0; iter < NUMITER; iter++) {

    /* 2.1. PROCESS ELEMENTS */
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int n_threads = omp_get_num_threads();
      int position = (tid+1)*(TOTALSIZE/n_threads);
      double saved_val;

      if(position < TOTALSIZE) {
        saved_val = V[position];
      }
	  #pragma omp barrier

      #pragma omp for
      for(i = 0; i < TOTALSIZE-1; i++) {
        if(i < position - 1) {
           V[i] = f(V[i], V[i+1]);
        } else {
           V[position - 1] = f(V[position - 1], saved_val); 
        }
      }

    }
    
    /* 2.2. END ITERATIONS LOOP */
  }

  /* 3. OUTPUT FINAL VALUES */
if(1) {
  printf("Output:\n"); 
  for(i = 0; i < TOTALSIZE; i++) {
    printf("%4d %f\n", i, V[i]);
  }
}

}
