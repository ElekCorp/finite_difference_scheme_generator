#!/bin/bash

#SBATCH --array=1-500
#SBATCH --job-name=TestJob
#SBATCH --output=TestJob.out
#SBATCH --time=1-00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M

####mkdir tmp

###  zig build-exe  main.zig -lc -O ReleaseFast
###  gcc -O3 -march=native -flto -Wall -o main_gcc spectral_bench.cpp -lm -flto
###  gcc -Ofast -march=native -flto -Wall -o main_gcc_ff spectral_bench.cpp -lm -flto

n=$(((($SLURM_ARRAY_TASK_ID + 1) * 1) % 50 ))
    start=$(date +%s.%N)
    python scheme_gen.py $n
    dur0=$(echo "$(date +%s.%N) - $start" | bc)

    start=$(date +%s.%N)
    ./main1 $n
    dur1=$(echo "$(date +%s.%N) - $start" | bc)

    start=$(date +%s.%N)
    ./main2 $n
    dur2=$(echo "$(date +%s.%N) - $start" | bc)

    ###printf "%s\n%.6f seconds\n%.6f seconds\n" $n $dur0 $dur1  > ./tmp/$((SLURM_ARRAY_TASK_ID)).txt
####result=`expr "$dur0 / $dur1" | bc -l`
####printf "%.6f\n" $result > ./tmp/$((SLURM_ARRAY_TASK_ID)).txt

printf "%i\n%.10f\n%.10f\n%.10f\n" $n $dur0 $dur1 $dur2 > ./tmp/$((SLURM_ARRAY_TASK_ID)).txt
