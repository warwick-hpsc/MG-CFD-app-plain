#!/bin/bash

#BSUB -J MG-CFD.<RUN ID>
#BSUB -o lsf.stdout
#BSUB -e lsf.stderr

#BSUB -nnodes 1
#BSUB -W <HOURS>:<MINUTES>
#BSUB -x
#BSUB -q <PARTITION>
