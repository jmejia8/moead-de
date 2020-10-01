CC=g++
CFLAGS=-I.

main : main.cpp
	$(CC)  main.cpp -o main
test :
	./main
clean :
	rm main
