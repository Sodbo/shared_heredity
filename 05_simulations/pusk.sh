#!/bin/bash

#SBATCH --job-name=testing
#SBATCH --output=_slurm_testing.out
#SBATCH --nodelist=mga-n3

srun R -f zapusk.R
