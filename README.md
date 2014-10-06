CPD
===
Testing the "cost" function computation time:

  1) Compile test.c with math libraby:
gcc test.c -o test -lm

  2) Run unix time to compare the times of computing or not the function:
time ./test t
time ./test

  The argument t makes the program execute the function.
