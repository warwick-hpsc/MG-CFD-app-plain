#!/bin/bash

#BSUB -J MG-CFD
#BSUB -o lsf.stdout
#BSUB -e lsf.stderr

#BSUB -nnodes 1
#BSUB -W <HOURS>:<MINUTES>
#BSUB -x
#BSUB -q <PARTITION>

cd <RUN_DIR>
