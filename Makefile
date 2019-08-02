#####################
### Macro options ###
###########################
## General settings
# LOG - print out logging messages

## Monitoring:
# TIME - monitor runtime of kernels
# PAPI - monitor PAPI event counters for each kernel

## Compute optimisation:
# OPT_LEVEL=x   - optimisation level. Can be 0, 1, 2 or 3
# OMP           - enable strong scaling with OpenMP threads
# OMP_SCATTERS  - allow multithreading of kernels containing scatter writes
# SIMD          - enable vectorisation (using OMP directives)
# DBLS_PER_SIMD - specify SIMD width
# FLUX_FISSION  - separate out the indirect writes from compute_flux_edge(), put into update()

## Data access
# LEGACY_ORDERING - Order edges to match output of Warwick's Python HDF5 conversion script. Used for internal validation


#################################
## Construct compilation command:
#################################

WARNINGS := -w

ifeq ($(COMPILER),gnu)
	CPP := g++
	CFLAGS += -fopenmp
	CFLAGS += -fmax-errors=1

	## Ensure GNU vectorizer is enabled
	OPTIMISATION += -ftree-vectorize

	GCC_OPT_REPORT_OPTIONS := 
	CFLAGS += $(GCC_OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -march=native
	CPU_SSE41_EXEC_TARGET = -msse4.1
	CPU_SSE42_EXEC_TARGET = -msse4.2
	CPU_AVX_EXEC_TARGET = -mavx
	CPU_AVX2_EXEC_TARGET = -mavx2
	KNL_AVX512_EXEC_TARGET = -mavx512f -mavx512er -mavx512cd -mavx512pf -march=knl
	CPU_AVX512_EXEC_TARGET = -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi -march=skylake-avx512

else ifeq ($(COMPILER),intel)
	CPP := icpc
	CFLAGS += -qopenmp
	CFLAGS += -fmax-errors=1

	INTEL_OPT_REPORT_OPTIONS := -qopt-report-phase=all -qopt-report=5
	CFLAGS += $(INTEL_OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -xHost
	CPU_SSE41_EXEC_TARGET = -xSSE4.1
	CPU_SSE42_EXEC_TARGET = -xSSE4.2
	CPU_AVX_EXEC_TARGET = -xAVX
	CPU_AVX2_EXEC_TARGET = -xCORE-AVX2
	KNL_AVX512_EXEC_TARGET = -xMIC-AVX512 -qopt-zmm-usage=high
	CPU_AVX512_EXEC_TARGET = -xCORE-AVX512 -qopt-zmm-usage=high

else ifeq ($(COMPILER),clang)
	CPP := clang++
	CFLAGS += -fopenmp
	CFLAGS += -fmax-errors=1

	OPT_REPORT_OPTIONS := -Rpass-missed=loop-vec -Rpass-analysis=loop-vec
	CFLAGS += $(OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -march=native
	CPU_SSE41_EXEC_TARGET = -msse4.1
	CPU_SSE42_EXEC_TARGET = -msse4.2
	CPU_AVX_EXEC_TARGET = -mavx
	CPU_AVX2_EXEC_TARGET = -mavx2
	KNL_AVX512_EXEC_TARGET = -mavx512f -mavx512er -mavx512cd -mavx512pf -march=knl
	CPU_AVX512_EXEC_TARGET = -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi -march=skylake-avx512

else ifeq ($(COMPILER),cray)
	CPP := CC

	WARNINGS := 

	OPT_REPORT_OPTIONS := -hlist=a
	CFLAGS += $(OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = 
	# Cray does not support Intel architectures older than Sandy Bridge, so cannot 
	# target SSE4.x 
	# CPU_SSE41_EXEC_TARGET = -target-cpu=barcelona
	# CPU_SSE42_EXEC_TARGET = -target-cpu=barcelona
	CPU_AVX_EXEC_TARGET = -target-cpu=sandybridge
	CPU_AVX2_EXEC_TARGET = -target-cpu=haswell
	CPU_AVX512_EXEC_TARGET = -target-cpu=skylake
	KNL_AVX512_EXEC_TARGET = -target-cpu=mic-knl

else
$(error Compiler not specified, aborting. Set 'COMPILER' to either "intel", "gnu", "clang" or "cray")
endif
CFLAGS += $(WARNINGS)

ifdef CPP_WRAPPER
	CPP := $(CPP_WRAPPER)
endif

ifdef OPT_LEVEL
	OPTIMISATION := -O$(OPT_LEVEL)
else
	OPTIMISATION := -O3
endif

## Disable aggressive floating-point optimization, to improve ability 
## of MG-CFD to assess floating-point performance. This is also necessary 
## to validate result against pre-computed solution.
ifeq ($(COMPILER),gnu)
	OPTIMISATION += -fno-fast-math
else ifeq ($(COMPILER),intel)
	OPTIMISATION += -fp-model precise
else ifeq ($(COMPILER),clang)
	OPTIMISATION += -fno-fast-math
else ifeq ($(COMPILER),cray)
	OPTIMISATION += -h fp2=noapprox
endif

ifdef INSN_SET
	ifeq ($(INSN_SET),Host)
		X_EXEC_CPU=$(HOST_EXEC_TARGET)
		X_EXEC_KNL=$(HOST_EXEC_TARGET)
	else ifeq ($(INSN_SET),SSE41)
		X_EXEC_CPU=$(CPU_SSE41_EXEC_TARGET)
		X_EXEC_KNL=$(CPU_SSE41_EXEC_TARGET)
	else ifeq ($(INSN_SET),SSE42)
		X_EXEC_CPU=$(CPU_SSE42_EXEC_TARGET)
		X_EXEC_KNL=$(CPU_SSE42_EXEC_TARGET)
	else ifeq ($(INSN_SET),AVX)
		X_EXEC_CPU=$(CPU_AVX_EXEC_TARGET)
		X_EXEC_KNL=$(CPU_AVX_EXEC_TARGET)
	else ifeq ($(INSN_SET),AVX2)
		X_EXEC_CPU=$(CPU_AVX2_EXEC_TARGET)
		X_EXEC_KNL=$(CPU_AVX2_EXEC_TARGET)
	else ifeq ($(INSN_SET),AVX512)
		X_EXEC_CPU=$(CPU_AVX512_EXEC_TARGET)
		X_EXEC_KNL=$(KNL_AVX512_EXEC_TARGET)
	else
		$(error Unknown value of 'INSN_SET')
	endif
else
	INSN_SET := Host
	X_EXEC_CPU=$(HOST_EXEC_TARGET)
	X_EXEC_KNL=$(HOST_EXEC_TARGET)
	BUILD_FLAGS += -DINSN_SET=$(INSN_SET)
endif

LIBS :=
ifneq (,$(findstring PAPI,$(BUILD_FLAGS)))
	LIBS += -lpapi -lpfm
endif

# BUILD_FLAGS_COMPRESSED := $(shell echo $(BUILD_FLAGS) | tr -d " ")
BUILD_FLAGS_COMPRESSED := $(COMPILER)$(shell echo $(BUILD_FLAGS) | tr -d " ")

# BIN_DIR = bin/$(shell hostname)
# OBJ_DIR = obj/$(shell hostname)/$(BUILD_FLAGS_COMPRESSED)
# OBJ_DIR_DBG = obj/$(shell hostname)/$(BUILD_FLAGS_COMPRESSED)_debug
BIN_DIR = bin
OBJ_DIR = obj/$(BUILD_FLAGS_COMPRESSED)
OBJ_DIR_DBG = obj/$(BUILD_FLAGS_COMPRESSED)_debug

#####################################################
## Now customise build to comply with BUILD_FLAGS: ##
#####################################################

SOURCES = src/euler3d_cpu_double.cpp \
		  src/Base/common.cpp \
		  src/Base/config.cpp \
		  src/Base/io.cpp \
		  src/Base/io_enhanced.cpp \
		  src/Kernels/flux_loops.cpp \
		  src/Kernels/cfd_loops.cpp \
		  src/Kernels/mg_loops.cpp \
		  src/Kernels/indirect_rw_loop.cpp \
		  src/Kernels/validation.cpp \
		  src/Monitoring/timer.cpp \
		  src/Monitoring/papi_funcs.cpp \
		  src/Monitoring/loop_stats.cpp \
		  src/Meshing/reorder.cpp \
		  src/Meshing/graph.cpp \
		  src/Meshing/colour.cpp \
		  src/Meshing/progress.cpp

OBJECTS     := $(patsubst src/%.cpp, $(OBJ_DIR)/%.o,     $(SOURCES))
OBJECTS_DBG := $(patsubst src/%.cpp, $(OBJ_DIR_DBG)/%.o, $(SOURCES))

INCLUDES := -Isrc -Isrc/Base -Isrc/Kernels -Isrc/Monitoring -Isrc/Meshing

#############
## TARGETS ##
#############

BIN_NAME := 
BIN_NAME_SUFFIX := $(BUILD_FLAGS_COMPRESSED)

main: BIN_NAME := $(BIN_DIR)/euler3d_cpu_double_$(BIN_NAME_SUFFIX).b
main: OPTIMISATION += $(X_EXEC_CPU)
main: generic

debug: BIN_NAME := $(BIN_DIR)/euler3d_cpu_double_debug_$(BIN_NAME_SUFFIX).b
debug: OPTIMISATION := -g -pg -O0
debug: BUILD_FLAGS += -DDEBUG
debug: generic_debug

$(OBJ_DIR)/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CPP) $(CFLAGS) $(OPTIMISATION) -c -o $@ $< $(BUILD_FLAGS) $(INCLUDES)

generic: $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(CPP) $(CFLAGS) $(OPTIMISATION) $^ -o $(BIN_NAME) $(BUILD_FLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR_DBG)/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CPP) $(CFLAGS) $(OPTIMISATION) -c -o $@ $< $(BUILD_FLAGS) $(INCLUDES)

generic_debug: $(OBJECTS_DBG)
	mkdir -p $(BIN_DIR)
	$(CPP) $(CFLAGS) $(OPTIMISATION) $^ -o $(BIN_NAME) $(BUILD_FLAGS) $(INCLUDES) $(LIBS)

clean:
	rm -rf $(BIN_DIR)
	# rm -rf obj/$(shell hostname)
	rm -rf obj
