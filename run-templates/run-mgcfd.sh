set -e
set -u

touch job-is-running.txt
if [ -f job-in-queue.txt ]; then
    rm job-in-queue.txt
fi

# Decide whether to compile app, execute app, or both:
do_execute=false
do_compile=false
if [ ! -z ${BATCH_EXECUTE_MODE+x} ]; then
  if [ "$BATCH_EXECUTE_MODE" = "execute" ]; then
    do_execute=true
  elif [ "$BATCH_EXECUTE_MODE" = "compile" ]; then
    do_compile=true
  fi
fi
for var in "$@" ; do
  if [ "$var" == "--execute" ]; then
    do_execute=true
  elif [ "$var" == "--compile" ]; then
    do_compile=true
  fi
done
if ! $do_execute ; then
  if ! $do_compile ; then
    # User must have executed this script without arguments. 
    # Assume that user expects script to both compile and execute app.
    do_compile=true
    do_execute=true
  fi
fi

# Compilation variables:
isa="<ISA>"
compiler="<COMPILER>"
cpp_wrapper="<CPP_WRAPPER>"
flags_final="<BUILD_FLAGS>"
debug=<DEBUG>

# File/dir paths:
run_outdir="<RUN_OUTDIR>"
parent_dir=`dirname "${run_outdir}"`
app_dirpath="<APP_DIRPATH>"
data_dirpath="<DATA_DIRPATH>"

# MG-CFD run variables:
_t=<NUM_THREADS>
_m=<MESH_MULTI>
mg_cycles=<MG_CYCLES>
validate_result=<VALIDATE_RESULT>


## Exit early if output csv files already exist.
if [ -f "${run_outdir}/Times.csv" ]; then
    echo "Times.csv already exists, meaning this job has already run."
    if [ -f "${run_outdir}"/job-is-running.txt ]; then
        rm "${run_outdir}"/job-is-running.txt
    fi
    exit 0
fi

# Compile:
if [ "$isa" != "" ]; then
    export INSN_SET="$isa"
    flags_final="${flags_final} -DINSN_SET=$isa"
fi
if $debug ; then
  bin_filename=euler3d_cpu_double_debug_"${compiler}"
else
  bin_filename=euler3d_cpu_double_"${compiler}"
fi
bin_filename="$bin_filename"`echo "$flags_final" | tr -d " "`.b
bin_filepath="${app_dirpath}/bin/${bin_filename}"

if $do_compile ; then
  export BUILD_FLAGS="$flags_final"
  cd "${app_dirpath}"
  make_cmd="COMPILER=${compiler} "
  if [ "$cpp_wrapper" != "" ]; then
    make_cmd+="CPP_WRAPPER=$cpp_wrapper "
  fi
  make_cmd+="make -j4 "
  if $debug ; then
    make_cmd+="debug"
  fi
  echo "$make_cmd"
  eval "$make_cmd"
  chmod a+x "$bin_filepath"
elif [ ! -f "$bin_filepath" ]; then
  echo "ERROR: Cannot find binary: $bin_filepath"
  rm "${run_outdir}"/job-is-running.txt
  exit 1
fi

# Grab object files:
if [ ! -d "${run_outdir}/objects" ]; then
    mkdir "${run_outdir}/objects"
fi
obj_dir="${app_dirpath}/obj/"
obj_dir+="${compiler}"
obj_dir+=`echo "$flags_final" | tr -d " "`
cp "${obj_dir}"/Kernels/flux_loops.o "${run_outdir}/objects/"
cp "${obj_dir}"/Kernels/indirect_rw_loop.o "${run_outdir}/objects/"


## Exit early if app execution not requested.
if ! $do_execute ; then
  echo "App execution not requested, exiting before app execution"
  exit 0
fi

# Execute:
if [ ! -z ${RUN_CMD+x} ]; then
  exec_command="$RUN_CMD "
else
  exec_command=""
fi
exec_command+="$bin_filepath -i input.dat -m $_m -p ${parent_dir}/papi.conf -o ${run_outdir}/ -g $mg_cycles"
if $validate_result ; then
  exec_command+=" -v"
fi
echo "EXECUTING $bin_filepath"
export OMP_NUM_THREADS=$_t
cd "${data_dirpath}"
echo "$exec_command"
eval "$exec_command"
rm "${run_outdir}"/job-is-running.txt
