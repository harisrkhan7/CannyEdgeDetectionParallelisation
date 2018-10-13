CC=clang++
all : move_upng main_o upng_o
	$(CC) --std=c++14 -o output.out main.o upng.o
main_o : main.cpp
	$(CC) -c main.cpp
upng_o: upng/upng.cpp
	$(CC) -c upng/upng.cpp
move_upng: 
	mv upng/upng.c upng/upng.cpp || true
run:
	./output.out
