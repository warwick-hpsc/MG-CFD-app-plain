#!/bin/bash

#MSUB -N MG-CFD
#MSUB -o moab.stdout
#MSUB -e moab.stderr

#MSUB -l nodes=1
#MSUB -l walltime=<HOURS>:<MINUTES>:00
#MSUB -n
#MSUB -q <PARTITION>

cd <RUN_DIR>

module load intel