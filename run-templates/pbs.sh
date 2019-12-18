#!/bin/bash

#PBS -N MG-CFD.<RUN ID>
#PBS -o pbs.stdout
#PBS -e pbs.stderr

#PBS -l select=1
#PBS -l walltime=<HOURS>:<MINUTES>:00
#PBS -A <BUDGET CODE>
#PBS -q <PARTITION>

#PBS -V

RUN_CMD="aprun -n 1 -d <NUM_THREADS>"

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
