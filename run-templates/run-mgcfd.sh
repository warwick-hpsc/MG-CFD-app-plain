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
app_dirpath="<APP_DIRPATH>"
data_dirpath="<DATA_DIRPATH>"

# MG-CFD run variables:
_t=<NUM_THREADS>
_m=<MESH_MULTI>
mg_cycles=<MG_CYCLES>
validate_result=<VALIDATE_RESULT>
renumber=<RENUMBER>
measure_mem_bound=<MEASURE_MEM_BOUND>
measure_compute_bound=<MEASURE_COMPUTE_BOUND>
run_synthetic_compute=<RUN_SYNTHETIC_COMPUTE>

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

objs_dirpath="${app_dirpath}"/obj/"${compiler}"`echo "$flags_final" | tr -d " "`

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

# Grab object files, compile logs, and optimisation reports:
if [ ! -d "${run_outdir}/objects" ]; then
    mkdir "${run_outdir}/objects"
fi
obj_dir="${app_dirpath}/obj/"
obj_dir+="${compiler}"
obj_dir+=`echo "$flags_final" | tr -d " "`
for kernel in flux unstructured_compute unstructured_stream compute_stream ; do
  for suffix in loop loops ; do
    # Compile logs:
    if [[ "$flags_final" = *"SIMD"* ]]; then
      log_fp="${obj_dir}/Kernels_vectorised/${kernel}_vec${suffix}".o.log
    else
      log_fp="${obj_dir}/Kernels/${kernel}_${suffix}".o.log
    fi
    echo "$log_fp"
    if [ -f "$log_fp" ]; then
      cp "$log_fp" "${run_outdir}"/objects/
    fi

    ## Optimisation reports:
    for ext in lst optrpt ; do
      if [[ "$flags_final" = *"SIMD"* ]]; then
        opt_fp="${obj_dir}/Kernels_vectorised/${kernel}_vec${suffix}.${ext}"
      else
        opt_fp="${obj_dir}/Kernels/${kernel}_${suffix}.${ext}"
      fi
      if [ -f "$opt_fp" ]; then
        cp "$opt_fp" "${run_outdir}"/objects/
      fi
    done

    # Objects:
    if [[ "$flags_final" = *"SIMD"* ]]; then
      obj_fp="${obj_dir}/Kernels_vectorised/${kernel}_vec${suffix}".o
    else
      obj_fp="${obj_dir}/Kernels/${kernel}_${suffix}".o
    fi
    if [ -f "$obj_fp" ]; then
      cp "$obj_fp" "${run_outdir}"/objects/
    fi

    # Update: run 'objdump' on the system to get assembly:
    if [[ "$flags_final" = *"SIMD"* ]]; then
      obj_fp="${run_outdir}"/objects/"${kernel}_vec${suffix}".o
    else
      obj_fp="${run_outdir}"/objects/"${kernel}_${suffix}".o
    fi
    if [ -f "$obj_fp" ]; then
      objdump_raw_command="objdump -d --no-show-raw-insn ${obj_fp}"
      objdump_raw_command+=" > ${obj_fp}.raw-asm"
      echo "$objdump_raw_command"
      eval "$objdump_raw_command"
      objdump_command="cat ${obj_fp}.raw-asm"
      objdump_command+=' | sed "s/^Disassembly of section/ fnc: Disassembly/g"'
      objdump_command+=' | sed "s/:$//g" | grep "^ " | grep ":" | sed "s/^[ \t]*//g"'
      objdump_command+=" > ${obj_fp}.asm"
      echo "$objdump_command"
      eval "$objdump_command"
    fi
  done
done

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
exec_command+="$bin_filepath -i input.dat -m $_m -p ${run_outdir}/papi.conf -o ${run_outdir}/ -g $mg_cycles"
if $validate_result ; then
  exec_command+=" -v"
fi
if $renumber ; then
  exec_command+=" -r"
fi
if $measure_mem_bound ; then
  exec_command+=" -b"
fi
if $measure_compute_bound ; then
  exec_command+=" -f"
fi
if $run_synthetic_compute ; then
  exec_command+=" -u"
fi

if [ "$compiler" = "intel" ]; then
  if [ -z ${KMP_AFFINITY+x} ]; then
    KMP_AFFINITY=""
  fi
  if [ -z ${KMP_GRANULARITY+x} ]; then
    KMP_GRANULARITY=""
  fi
  if [ "$KMP_AFFINITY" = "" ] && [ "$KMP_GRANULARITY" = "" ] ; then
    export KMP_AFFINITY=scatter
    export KMP_GRANULARITY=core
  fi
else
  if [ -z ${OMP_PROC_PLACES+x} ]; then
    OMP_PROC_PLACES=""
  fi
  if [ -z ${OMP_PROC_BIND+x} ]; then
    OMP_PROC_BIND=""
  fi
  if [ "$OMP_PROC_PLACES" = "" ] && [ "$OMP_PROC_BIND" = "" ] ; then
    export OMP_PLACES=sockets
    export OMP_PROC_BIND=spread
  fi
fi
if $renumber ; then
  exec_command+=" -r"
fi
echo "EXECUTING $bin_filepath"
export OMP_NUM_THREADS=$_t
cd "${data_dirpath}"
echo ""
echo "$exec_command"
eval "$exec_command"

touch "${run_outdir}"/job-is-complete.txt
rm "${run_outdir}"/job-is-running.txt
