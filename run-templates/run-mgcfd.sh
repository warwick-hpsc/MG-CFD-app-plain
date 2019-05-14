set -e
set -o

isa=<ISA>
run_outdir="<RUN_OUTDIR>"
parent_dir=`dirname "${run_outdir}"`
flags_final="<BUILD_FLAGS>"
app_dirpath="<APP_DIRPATH>"
_cc=<COMPILER>
_t=<NUM_THREADS>
data_dirpath="<DATA_DIRPATH>"
_m=<MESH_MULTI>
mg_cycles=<MG_CYCLES>

if [ -f "${run_outdir}/Times.csv" ]; then
    echo "Times.csv already exists, meaning this job has already run."
    exit 0
fi

# Compile:
if [ "$isa" != "" ]; then
    export INSN_SET="$isa"
    flags_final="${flags_final} -DINSN_SET=$isa"
fi
export BUILD_FLAGS="$flags_final"
cd "${app_dirpath}"
CC="$_cc" make -j4
bin_filename=euler3d_cpu_double_"${_cc}"
bin_filename="$bin_filename"`echo "$flags_final" | tr -d " "`.b

# Grab object files:
if [ ! -d "${run_outdir}/objects" ]; then
    mkdir "${run_outdir}/objects"
fi
obj_dir="${app_dirpath}/obj/`hostname`/"
obj_dir="${obj_dir}$_cc"
obj_dir="${obj_dir}"`echo "$flags_final" | tr -d " "`
cp "${obj_dir}"/Kernels/flux_kernels.o "${run_outdir}/objects/"
cp "${obj_dir}"/Kernels/indirect_rw_kernel.o "${run_outdir}/objects/"

# Execute:
cd "${data_dirpath}"
bin_filepath="${app_dirpath}/bin/`hostname`/${bin_filename}"
export OMP_NUM_THREADS=$_t
echo "EXECUTING $bin_filepath"
eval "$bin_filepath" -i input.dat -m $_m -p "${parent_dir}"/papi.conf -o "${run_outdir}/" -g $mg_cycles
