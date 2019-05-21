#!/bin/bash

#PBS -N MG-CFD.<RUN ID>
#PBS -o pbs.stdout
#PBS -e pbs.stderr

#PBS -l select=1
#PBS -l walltime=<HOURS>:<MINUTES>:00
#PBS -A <BUDGET CODE>
#PBS -q <PARTITION>

## Load compiler module(s):
compiler="<COMPILER>"
if [ "$compiler" = "intel" ]; then
  module swap PrgEnv-cray PrgEnv-intel
elif [ "$compiler" = "gnu" ]; then
  module swap PrgEnv-cray PrgEnv-gnu
# elif [ "$compiler" = "cray" ]; then
# elif [ "$compiler" = "clang" ]; then
fi

module load papi

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
