all:
	gcc test.c -o test -lm -fopenmp
clean:
	rm *~ *.o test

