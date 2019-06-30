set -e
set -u

touch job-is-running.txt
if [ -f job-in-queue.txt ]; then
    rm job-in-queue.txt
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
    rm "${run_outdir}"/job-is-running.txt
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
# bin_filepath="${app_dirpath}/bin/`hostname`/${bin_filename}"
bin_filepath="${app_dirpath}/bin/${bin_filename}"

# if [ ! -f "$bin_filepath" ]; then
    ## Try compiling anyway, source files may have changed
  if [[ `hostname` == *"login"* ]]; then
    ## On login node, compile
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
# fi

# Grab object files:
if [ ! -d "${run_outdir}/objects" ]; then
    mkdir "${run_outdir}/objects"
fi
# obj_dir="${app_dirpath}/obj/`hostname`/"
obj_dir="${app_dirpath}/obj/"
obj_dir+="${compiler}"
obj_dir+=`echo "$flags_final" | tr -d " "`
cp "${obj_dir}"/Kernels/flux_loops.o "${run_outdir}/objects/"
cp "${obj_dir}"/Kernels/indirect_rw_loop.o "${run_outdir}/objects/"

if [[ `hostname` == *"login"* ]]; then
  ## Assume on a login node, do not execute the code.
  echo "Detected presence on login node, aborting before app execution."
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
