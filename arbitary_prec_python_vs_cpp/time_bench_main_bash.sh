#!/bin/bash

#SBATCH --job-name=MainJob
#SBATCH --output=MainJob.out
#SBATCH --time=1-00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=500M

compID=$(sbatch --parsable time_test_compile.sh )
benchID=$(sbatch --parsable --dependency=afterany:${compID} time_test.sh)
PlotID=$(sbatch --parsable --dependency=afterany:${benchID} time_test_plot.sh)