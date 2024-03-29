all: clean
	gcc test.c -o test -lm -fopenmp
	g++ sequential.cpp -o seq -lm -fopenmp
	g++ parallel_first_try.cpp -o par -lm -fopenmp
clean:
	rm -f *~ *.o test main
run: all
	chmod +x tests.sh
	./tests.sh
