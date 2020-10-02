CC=g++
CFLAGS=-O0 -Wall
CFLAGSRELEASE=-O3

all : main moead.so

main : main.cpp
	$(CC)  main.cpp $(CFLAGS) -o main

moead.so : moead-de-lib.cpp
		$(CC)  moead-de-lib.cpp $(CFLAGS) -o moead-de-lib.so

test :
	./main

clean :
	rm main moead-de-lib.so
