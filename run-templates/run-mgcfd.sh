set -e
set -u

# Compilation variables:
isa="<ISA>"
compiler="<COMPILER>"
cpp_override="<CPP_OVERRIDE>"
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


if [ -f "${run_outdir}/Times.csv" ]; then
    echo "Times.csv already exists, meaning this job has already run."
    exit 0
fi


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
bin_filepath="${app_dirpath}/bin/`hostname`/${bin_filename}"

if [ ! -f "$bin_filepath" ]; then
	# Compile:
	export BUILD_FLAGS="$flags_final"
	cd "${app_dirpath}"
	make_cmd="COMPILER=${compiler} "
	if [ "$cpp_override" != "" ]; then
		make_cmd+="CPP_OVERRIDE=$cpp_override "
	fi
	make_cmd+="make -j4 "
	if $debug ; then
		make_cmd+="debug"
	fi
	echo "$make_cmd"
	eval "$make_cmd"
fi

if [[ `hostname` == *"login"* ]]; then
	## Assume on a login node, do not execute the code.
	echo "Detected presence on login node, aborting before app execution."
	exit 0
fi

# Grab object files:
if [ ! -d "${run_outdir}/objects" ]; then
    mkdir "${run_outdir}/objects"
fi
obj_dir="${app_dirpath}/obj/`hostname`/"
obj_dir="${obj_dir}${compiler}"
obj_dir="${obj_dir}"`echo "$flags_final" | tr -d " "`
cp "${obj_dir}"/Kernels/flux_loops.o "${run_outdir}/objects/"
cp "${obj_dir}"/Kernels/indirect_rw_loop.o "${run_outdir}/objects/"

# Execute:
cd "${data_dirpath}"
export OMP_NUM_THREADS=$_t
echo "EXECUTING $bin_filepath"
if $debug ; then
	gdb --args "$bin_filepath" -i input.dat -m $_m -p "${parent_dir}"/papi.conf -o "${run_outdir}/" -g $mg_cycles
elif $validate_result ; then
	eval "$bin_filepath" -i input.dat -m $_m -p "${parent_dir}"/papi.conf -o "${run_outdir}/" -g $mg_cycles -v
else
	eval "$bin_filepath" -i input.dat -m $_m -p "${parent_dir}"/papi.conf -o "${run_outdir}/" -g $mg_cycles
fi
