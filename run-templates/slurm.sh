#!/bin/bash

#SBATCH --job-name=MG-CFD.<RUN ID>
#SBATCH --output=sbatch.stdout
#SBATCH --error=sbatch.stderr

#SBATCH --nodes=1
#SBATCH --time=<HOURS>:<MINUTES>:00
#SBATCH --exclusive
#SBATCH --partition=<PARTITION>
#SBATCH --account=<BUDGET CODE>

#SBATCH --export=ALL

RUN_CMD="srun"
