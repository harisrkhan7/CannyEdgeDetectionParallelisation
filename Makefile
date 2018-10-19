CC=mpicc
all : move_upng main_o upng_o
	$(CC) --std=c++14 -lstdc++ -lomp -o output.out main.o upng.o
main_o : main.cpp
	$(CC) -Xpreprocessor -fopenmp -c main.cpp
upng_o: upng/upng.cpp
	$(CC) -c upng/upng.cpp
move_upng: 
	mv upng/upng.c upng/upng.cpp || true
run:
	./output.out img/lion.png
mpirun:
	OMPI_MCA_btl=^vader mpirun -np 4 ./output.out img/lion.png
