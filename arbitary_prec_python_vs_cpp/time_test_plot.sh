#!/bin/bash

#SBATCH --job-name=PlotJob
#SBATCH --output=PlotJob.out
#SBATCH --time=1-00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=500M

python time_bench_plot.py 500

rm -rf main main1 main2
