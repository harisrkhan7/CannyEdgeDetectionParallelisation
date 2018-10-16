#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --time=0-12:00:00

module load OpenMPI/3.0.0-GCC-6.2.0
module load GCC/6.2.0

make

num_threads_arr=("1" "2" "3" "4")
filename_arr=("img/smallwheel.png" "img/wheel.png" "img/baboon.png" "img/elephant.png" "img/lion.png" "img/gazelle.png" "img/giraffe.png")

echo "filename,width,height,size,omp_num_threads,time(s)"
for img_name in "${filename_arr[@]}"; do
    for num_threads in "${num_threads_arr[@]}"; do
        elapsed_time="$(OMP_NUM_THREADS=$num_threads ./output.out $img_name 2> /dev/null)"

        echo -n "$img_name"
        echo -n ","
        img_width="$(head -n 3 out/out.pgm | tail -1)"
        echo -n $img_width
        echo -n ","
        img_height="$(head -n 2 out/out.pgm | tail -1)"
        echo -n $img_height
        echo -n ","
        echo -n $num_threads
        echo -n ","
        echo $elapsed_time
        
    done
done
