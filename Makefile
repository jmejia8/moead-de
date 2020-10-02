CC=g++
CFLAGS=-O3

# main : main.cpp
# 	$(CC)  main.cpp $(CFLAGS) -o main

moead.so : moead-de-lib.cpp
		$(CC)  moead-de-lib.cpp $(CFLAGS) -o moead-de-lib.so
test :
	./main
clean :
	rm main
