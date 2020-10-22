#!/bin/bash

#MSUB -N MG-CFD.<RUN ID>
#MSUB -o moab.stdout
#MSUB -e moab.stderr

#MSUB -l nodes=1
#MSUB -l procs=1
#MSUB -l walltime=<HOURS>:<MINUTES>:00
#MSUB -n
#MSUB -q <PARTITION>
#MSUB -A <BUDGET CODE>

#MSUB -V

export RUN_CMD="srun --cpus-per-task=<NUM_THREADS>"
