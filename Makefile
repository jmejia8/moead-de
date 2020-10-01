CC=g++
CFLAGS=-O3

main : main.cpp
	$(CC)  main.cpp $(CFLAGS) -o main
test :
	./main
clean :
	rm main
