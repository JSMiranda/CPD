all: clean
	gcc test.c -o test -lm -fopenmp
	g++ sequential.cpp -o main -lm -fopenmp
clean:
	rm -f *~ *.o test main
run: all
	chmod +x tests.sh
	./tests.sh
