CC=clang++
all : main_o upng_o
	cc --std=c++14 -o a.out main.o upng.o
main_o : main.cpp 
	cc -c main.cpp
upng_o: upng/upng.cpp
	cc -c upng/upng.cpp

