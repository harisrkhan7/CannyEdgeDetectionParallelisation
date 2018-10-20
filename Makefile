CC=mpic++
all : move_upng main_o upng_o
	$(CC) --std=c++14 -lstdc++ -lgomp -lm -o output.out main.o upng.o
main_o : main.cpp
	$(CC) -Xpreprocessor -fopenmp -c main.cpp
upng_o: upng/upng.cpp
	$(CC) -c upng/upng.cpp
move_upng: 
	mv upng/upng.c upng/upng.cpp || true && mkdir out || true
run:
	./output.out img/lion.png
mpirun:
	mpirun -np 4 ./output.out img/smallwheel.png
