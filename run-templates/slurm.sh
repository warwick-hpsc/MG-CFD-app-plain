#!/bin/bash

#SBATCH --job-name=MG-CFD
#SBATCH --output=sbatch.stdout
#SBATCH --error=sbatch.stderr

#SBATCH --nodes=1
#SBATCH --time=<HOURS>:<MINUTES>:00
#SBATCH --exclusive
#SBATCH --partition=<PARTITION>

cd <RUN_DIR>

module load intel