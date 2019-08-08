#!/bin/bash

#BSUB -J MG-CFD.<RUN ID>
#BSUB -o lsf.stdout
#BSUB -e lsf.stderr

#BSUB -nnodes 1
#BSUB -R "affinity[core(<NUM_THREADS>, same=socket)]"
#BSUB -W <HOURS>:<MINUTES>
#BSUB -x
#BSUB -q <PARTITION>
#BSUB -P <BUDGET CODE>

## Note: LSF documentation claims that user environment variables 
##       are exported by default. However I have not confirmed this.
##       This feature is necessary, so fingers crossed!

RUN_CMD="mpirun -n 1"
