#!/bin/bash

#PBS -N MG-CFD
#PBS -o pbs.stdout
#PBS -e pbs.stderr

#PBS -l select=1
#PBS -l walltime=<HOURS>:<MINUTES>:00
#PBS -A <BUDGET CODE>
#PBS -q <PARTITION>

cd <RUN_DIR>

module load papi
echo "CRAY_LD_LIBRARY_PATH = $CRAY_LD_LIBRARY_PATH"
if [ "$CRAY_LD_LIBRARY_PATH" != "" ]; then
	if [ "$LD_LIBRARY_PATH" != "" ]; then
		export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${CRAY_LD_LIBRARY_PATH}"
	else
		export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}"
	fi

	if [ "$LIBRARY_PATH" != "" ]; then
	  export LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
	else
	  export LIBRARY_PATH="$LD_LIBRARY_PATH"
	fi
fi
if [ "$CPATH" != "" ]; then
	export CPATH="${CPATH}:/opt/cray/papi/5.5.1.1/include"
else
	export CPATH="/opt/cray/papi/5.5.1.1/include"
fi

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
