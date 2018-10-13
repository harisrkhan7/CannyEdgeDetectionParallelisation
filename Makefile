CC=clang++
all : main_o upng_o
	$(CC) --std=c++14 -o a.out main.o upng.o
main_o : main.cpp 
	$(CC) -c main.cpp
upng_o: upng/upng.cpp
	$(CC) -c upng/upng.cpp
