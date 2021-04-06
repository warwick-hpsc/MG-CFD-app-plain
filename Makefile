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

ifeq ($(COMPILER),gnu)
	CPP = g++
else ifeq ($(COMPILER),intel)
	CPP = icpc
else ifeq ($(COMPILER),clang)
	CPP = clang++
else ifeq ($(COMPILER),cray)
	CPP = CC
else ifeq ($(COMPILER),arm)
	CPP = armclang++
else ifeq ($(COMPILER),fujitsu)
	CPP = FCC
else
$(error Compiler not specified, aborting. Set 'COMPILER' to either "intel", "gnu", "clang" or "cray")
endif
COMPILER_API := $(COMPILER)

ifeq ($(COMPILER),cray)
	## Check whether Cray uses Clang frontend:
	_v = $(shell CC --help 2>/dev/null | grep Clang | head -n1 | grep -o Clang)
	ifeq ($(_v),Clang)
	    # Yes, this Cray does just wrap Clang.
	    COMPILER_API = clang
	endif
else ifeq ($(COMPILER),arm)
	COMPILER_API = clang
endif

## Handle code optimisation:
ifdef OPT_LEVEL
	OPTIMISATION := -O$(OPT_LEVEL)
else
	OPTIMISATION := -O3
endif
ifeq (DPRECISE_FP,$(findstring DPRECISE_FP, $(BUILD_FLAGS)))
	PRECISE_FP = yes
else
	PRECISE_FP = no
endif
ifeq ($(PRECISE_FP),yes)
	ifeq ($(COMPILER_API),gnu)
		OPTIMISATION += -fno-fast-math
				
		## Disable C math function error checking, as prevents SIMD:
		OPTIMISATION += -fno-math-errno
		
	else ifeq ($(COMPILER_API),intel)
		OPTIMISATION += -fp-model precise
		
	else ifeq ($(COMPILER_API),clang)
		OPTIMISATION += -fno-fast-math

		## Disable C math function error checking, as prevents SIMD:
		OPTIMISATION += -fno-math-errno
		
	else ifeq ($(COMPILER_API),cray)
		OPTIMISATION += -h fp2=noapprox
		
	else ifeq ($(COMPILER_API),fujitsu)
		OPTIMISATION += -Knofp_relaxed
	endif
else
	ifeq ($(COMPILER_API),gnu)
		OPTIMISATION += -ffast-math

	else ifeq ($(COMPILER_API),intel)
		OPTIMISATION += -fp-model fast=2

	else ifeq ($(COMPILER_API),clang)
		OPTIMISATION += -ffast-math

	else ifeq ($(COMPILER_API),cray)
		OPTIMISATION += -h fp2=approx

	else ifeq ($(COMPILER_API),fujitsu)
		OPTIMISATION += -Kfp_relaxed
	endif
endif
ifeq (DSIMD,$(findstring DSIMD, $(BUILD_FLAGS)))
	SIMD = yes
else
	SIMD = no
endif
ifeq ($(SIMD),yes)
	ifeq ($(COMPILER_API),gnu)
		OPTIMISATION += -ftree-vectorize

	else ifeq ($(COMPILER_API),intel)
		
	else ifeq ($(COMPILER_API),clang)
		## TODO: review this, is it still needed?
		# OPTIMISATION += -mllvm -vectorizer-min-trip-count=2

	else ifeq ($(COMPILER_API),cray)

	else ifeq ($(COMPILER_API),fujitsu)
		OPTIMISATION += -Ksimd 
		OPTIMISATION += -Ksimd_packed_promotion
		# OPTIMISATION += -fopenmp-simd
	endif
else
	ifeq ($(COMPILER_API),gnu)

	else ifeq ($(COMPILER_API),intel)

	else ifeq ($(COMPILER_API),clang)

	else ifeq ($(COMPILER_API),cray)

	else ifeq ($(COMPILER_API),fujitsu)
		OPTIMISATION += -Knosimd -KNOSVE
	endif
endif

WARNINGS := -w

ifeq ($(COMPILER_API),gnu)
	CFLAGS += -fopenmp
	CFLAGS += -fmax-errors=1

	# ## Enable all warnings, and treat as errors, to help cleanup code:
	# WARNINGS := -Wall -Wunused-function -Wunused-parameter -Werror

	HOST_EXEC_TARGET = -march=native
	CPU_SSE41_EXEC_TARGET = -msse4.1
	CPU_SSE42_EXEC_TARGET = -msse4.2
	CPU_AVX_EXEC_TARGET = -mavx
	CPU_AVX2_EXEC_TARGET = -mavx2
	KNL_AVX512_EXEC_TARGET = -mavx512f -mavx512er -mavx512cd -mavx512pf -march=knl
	CPU_AVX512_EXEC_TARGET = -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi -march=skylake-avx512
	AARCH64_EXEC_TARGET = -march=armv8-a
	A64FX_EXEC_TARGET = -mcpu=a64fx

else ifeq ($(COMPILER_API),intel)
	CFLAGS += -qopenmp
	CFLAGS += -fmax-errors=1
	CFLAGS += -vec-threshold0

	INTEL_OPT_REPORT_OPTIONS := -qopt-report-phase=all -qopt-report=5
	# INTEL_OPT_REPORT_OPTIONS := -qopt-report-phase=vec -qopt-report=5
	CFLAGS += $(INTEL_OPT_REPORT_OPTIONS)

	# ## Enable all warnings, and treat as errors, to help cleanup code:
	# WARNINGS := -Wall -Wunused-function -Wunused-parameter -Werror
	# ## ... except for when I explicitly disable SIMD:
	# 	# Diagnostic code 15552 is: error #15552: loop was not vectorized with "simd"
	# 	# This warning can be ignored, as it only appears when simd is explicitly blocked 
	# 	# with a pragma: #pragma omp simd safelen(1)
	# WARNINGS += -diag-disable=15552

	HOST_EXEC_TARGET = -xHost
	CPU_SSE41_EXEC_TARGET = -xSSE4.1
	CPU_SSE42_EXEC_TARGET = -xSSE4.2
	CPU_AVX_EXEC_TARGET = -xAVX
	CPU_AVX2_EXEC_TARGET = -xCORE-AVX2
	KNL_AVX512_EXEC_TARGET = -xMIC-AVX512 -qopt-zmm-usage=high
	CPU_AVX512_EXEC_TARGET = -xCORE-AVX512 -qopt-zmm-usage=high

else ifeq ($(COMPILER_API),clang)
	CFLAGS += -fopenmp
	CFLAGS += -ferror-limit=1
	CFLAGS += -finline-hint-functions

	## Loop unroller interferes with vectorizer, disable:
	OPTIMISATION += -fno-unroll-loops

	OPT_REPORT_OPTIONS := 
	OPT_REPORT_OPTIONS += -Rpass-missed=loop-vec ## Report SIMD failures
	OPT_REPORT_OPTIONS += -Rpass="loop-(unroll|vec)" ## Report loop transformations
	OPT_REPORT_OPTIONS += -Rpass-analysis=loop-vectorize ## Report WHY vectorize failed
	OPT_REPORT_OPTIONS += -fsave-optimization-record -gline-tables-only -gcolumn-info
	CFLAGS += $(OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -mcpu=native
	CPU_SSE41_EXEC_TARGET = -msse4.1
	CPU_SSE42_EXEC_TARGET = -msse4.2
	CPU_AVX_EXEC_TARGET = -mavx
	CPU_AVX2_EXEC_TARGET = -mavx2
	KNL_AVX512_EXEC_TARGET = -mavx512f -mavx512er -mavx512cd -mavx512pf -march=knl
	CPU_AVX512_EXEC_TARGET = -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi -march=skylake-avx512
	AARCH64_EXEC_TARGET = --target=aarch64
	A64FX_EXEC_TARGET = -mcpu=a64fx

else ifeq ($(COMPILER_API),cray)
	## Loop unroller interferes with vectorizer, disable:
	## Todo: how does non-LLVM Cray expose loop unroller?

    CFLAGS += -h omp

	WARNINGS := 

	OPT_REPORT_OPTIONS :=
	OPT_REPORT_OPTIONS += -hlist=a
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
	AARCH64_EXEC_TARGET = 

else ifeq ($(COMPILER_API),fujitsu)
	## Until I figure out what predefined macros exist, create my own:
	CFLAGS += -D__COMPILER_FCC__

	OPTIMISATION += -Kfast

	OPT_REPORT_OPTIONS = -Nlst=t
	CFLAGS += $(OPT_REPORT_OPTIONS)

	A64FX_EXEC_TARGET = -KA64FX 

else
$(error Compiler not specified, aborting. Set 'COMPILER' to either "intel", "gnu", "clang" or "cray")
endif
CFLAGS += $(WARNINGS)

ifdef CPP_WRAPPER
	CPP := $(CPP_WRAPPER)
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
else ifeq ($(INSN_SET),AARCH64)
	X_EXEC_CPU=$(AARCH64_EXEC_TARGET)
else ifeq ($(INSN_SET),A64FX)
	X_EXEC_CPU=$(A64FX_EXEC_TARGET)
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
	ifeq ($(COMPILER),gnu)
        ifdef PAPI_INCLUDE_PATH
    		INCLUDES += -I$(PAPI_INCLUDE_PATH)
        endif
	else ifeq ($(COMPILER),cray)
        ifdef PAPI_LIB_PATH
    		LIBS += -L$(PAPI_LIB_PATH)
        endif
	endif
	LIBS += -lpapi -lpfm
endif

BUILD_FLAGS_COMPRESSED := $(COMPILER)$(shell echo $(BUILD_FLAGS) | tr -d " ")

BIN_DIR = bin
OBJ_DIR = obj/$(BUILD_FLAGS_COMPRESSED)
OBJ_DIR_DBG = obj/$(BUILD_FLAGS_COMPRESSED)_debug

ifeq ($(COMPILER_API),gnu)
	OPT_REPORT_OPTIONS = 
	OPT_REPORT_OPTIONS += -fopt-info-vec-all=$(OBJ_DIR)/gcc-opt-info
	CFLAGS += $(OPT_REPORT_OPTIONS)
endif

#####################################################
## Now customise build to comply with BUILD_FLAGS: ##
#####################################################

SOURCES = src/euler3d_cpu_double.cpp \
		  src/Base/common.cpp \
		  src/Base/config.cpp \
		  src/Base/io.cpp \
		  src/Base/io_enhanced.cpp \
		  src/kernel_wrappers.cpp \
		  src/Kernels/flux_loops.cpp \
		  src/Kernels/cfd_loops.cpp \
		  src/Kernels/mg_loops.cpp \
		  src/Kernels/unstructured_stream_loop.cpp \
		  src/Kernels/unstructured_compute_loop.cpp \
		  src/Kernels/compute_stream_loop.cpp \
		  src/Kernels_vectorised/flux_vecloops.cpp \
		  src/Kernels_vectorised/unstructured_stream_vecloop.cpp \
		  src/Kernels_vectorised/unstructured_compute_vecloop.cpp \
		  src/Kernels_vectorised/compute_stream_vecloop.cpp \
		  src/Kernels/validation.cpp \
		  src/Monitoring/timer.cpp \
		  src/Monitoring/papi_funcs.cpp \
		  src/Monitoring/loop_stats.cpp \
		  src/Meshing/reorder.cpp \
		  src/Meshing/graph.cpp \
		  src/Meshing/colour.cpp \
		  src/Meshing/reduce_bw.cpp \
		  src/Meshing/progress.cpp

OBJECTS     := $(patsubst src/%.cpp, $(OBJ_DIR)/%.o,     $(SOURCES))
OBJECTS_DBG := $(patsubst src/%.cpp, $(OBJ_DIR_DBG)/%.o, $(SOURCES))

INCLUDES += -Isrc -Isrc/Base -Isrc/Kernels -Isrc/Kernels_vectorised -Isrc/Monitoring -Isrc/Meshing

#############
## TARGETS ##
#############

BIN_NAME := 
BIN_NAME_SUFFIX := $(BUILD_FLAGS_COMPRESSED)

main: BIN_NAME := $(BIN_DIR)/euler3d_cpu_double_$(BIN_NAME_SUFFIX).b
main: OPTIMISATION += $(X_EXEC_CPU)
main: clean_reports
main: generic

debug: BIN_NAME := $(BIN_DIR)/euler3d_cpu_double_debug_$(BIN_NAME_SUFFIX).b
debug: OPTIMISATION := -g -pg -O0
debug: BUILD_FLAGS += -DDEBUG
debug: generic_debug

clean_reports:
	rm -f *.lst
	rm -f gcc-report*

$(OBJ_DIR)/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CPP) $(CFLAGS) $(OPTIMISATION) -c -o $@ $< $(BUILD_FLAGS) $(INCLUDES) 2>&1 | tee $@.log
# $(CPP) $(CFLAGS) $(OPTIMISATION) -c -o $@ $< $(BUILD_FLAGS) $(INCLUDES) > $@.log 2>&1

generic: $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(CPP) $(CFLAGS) $(OPTIMISATION) $^ -o $(BIN_NAME) $(BUILD_FLAGS) $(INCLUDES) $(LIBS)
	if ls *.lst 2>/dev/null ; then mv *.lst $(OBJ_DIR)/ ; fi

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
