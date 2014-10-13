all:
	gcc test.c -o test -lm -fopenmp
	g++ sequential.cpp -o main -lm -fopenmp
clean:
	rm -f *~ *.o test main
run:
	./main public-instances/ex10.15.in
