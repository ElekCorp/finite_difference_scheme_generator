#!/bin/bash

#SBATCH --job-name=compile
#SBATCH --output=compile.out
#SBATCH --time=1-00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M

if [ ! -d "./tmp"  ]
then 
    mkdir tmp
    
fi


###zig build-exe  main.zig -lc -O ReleaseFast
###zig cc -O3 -Wall -o main main.cpp -lgmp -lgmpxx
###zig cc -Ofast -march=native -flto -Wall -o main1 main.cpp -lm -flto -lgmp -lgmpxx
g++ -O0 -march=native -flto -Wall -o main1 main.cpp -lm -flto -lgmp -lgmpxx -lstdc++
g++ -O3 -march=native -flto -Wall -o main2 main.cpp -lm -flto -lgmp -lgmpxx -lstdc++
