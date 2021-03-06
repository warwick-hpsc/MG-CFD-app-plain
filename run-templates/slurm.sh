#!/bin/bash

#SBATCH --job-name=MG-CFD.<RUN ID>
#SBATCH --output=sbatch.stdout
#SBATCH --error=sbatch.stderr

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<NUM_THREADS>
#SBATCH --time=<HOURS>:<MINUTES>:00
#SBATCH --exclusive
#SBATCH --partition=<PARTITION>
#SBATCH --account=<BUDGET CODE>

#SBATCH --export=ALL

export RUN_CMD="srun --cpu_bind=cores"
