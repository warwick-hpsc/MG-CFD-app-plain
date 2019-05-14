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

ifeq ($(CC),gnu)
	CPP := g++
	CPP += -fopenmp

	GCC_OPT_REPORT_OPTIONS := 
	CPP += $(GCC_OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -march=native
	CPU_SSE41_EXEC_TARGET = -msse4.1
	CPU_SSE42_EXEC_TARGET = -msse4.2
	CPU_AVX_EXEC_TARGET = -mavx
	CPU_AVX2_EXEC_TARGET = -mavx2
	KNL_AVX512_EXEC_TARGET = -mavx512f -mavx512er -mavx512cd -mavx512pf -march=knl
	CPU_AVX512_EXEC_TARGET = -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi -march=skylake-avx512

else ifeq ($(CC),intel)
	CPP := icpc
	CPP += -qopenmp

	INTEL_OPT_REPORT_OPTIONS := 
	CPP += $(INTEL_OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -xHost
	CPU_SSE41_EXEC_TARGET = -xSSE4.1
	CPU_SSE42_EXEC_TARGET = -xSSE4.2
	CPU_AVX_EXEC_TARGET = -xAVX
	CPU_AVX2_EXEC_TARGET = -xCORE-AVX2
	KNL_AVX512_EXEC_TARGET = -xMIC-AVX512 -qopt-zmm-usage=high
	CPU_AVX512_EXEC_TARGET = -xCORE-AVX512 -qopt-zmm-usage=high

else ifeq ($(CC),clang)
	CPP := clang++
	CPP += -fopenmp

	OPT_REPORT_OPTIONS := 
	CPP += $(OPT_REPORT_OPTIONS)

	HOST_EXEC_TARGET = -march=native
	CPU_SSE41_EXEC_TARGET = -msse4.1
	CPU_SSE42_EXEC_TARGET = -msse4.2
	CPU_AVX_EXEC_TARGET = -mavx
	CPU_AVX2_EXEC_TARGET = -mavx2
	KNL_AVX512_EXEC_TARGET = -mavx512f -mavx512er -mavx512cd -mavx512pf -march=knl
	CPU_AVX512_EXEC_TARGET = -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi -march=skylake-avx512

else
$(error Compiler not specified, aborting. Set 'CC' to either "intel", "gnu" or "clang")
endif
CPP += $(WARNINGS)
CPP += -fmax-errors=1

ifdef OPT_LEVEL
	OPTIMISATION := -O$(OPT_LEVEL)
else
	OPTIMISATION := -O3
endif

## Disable aggressive floating-point optimization, to improve ability 
## of MG-CFD to assess floating-point performance
ifeq ($(CC),gnu)
	OPTIMISATION += -fno-fast-math
else ifeq ($(CC),intel)
	OPTIMISATION += -fp-model precise
else ifeq ($(CC),clang)
	OPTIMISATION += -fno-fast-math
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
endif

LIBS :=
ifneq (,$(findstring PAPI,$(BUILD_FLAGS)))
	LIBS += -lpapi -lpfm
endif

# BUILD_FLAGS_COMPRESSED := $(shell echo $(BUILD_FLAGS) | tr -d " ")
BUILD_FLAGS_COMPRESSED := $(CC)$(shell echo $(BUILD_FLAGS) | tr -d " ")

BIN_DIR = bin/$(shell hostname)
OBJ_DIR = obj/$(shell hostname)/$(BUILD_FLAGS_COMPRESSED)
OBJ_DIR_DBG = obj/$(shell hostname)/$(BUILD_FLAGS_COMPRESSED)_debug

#####################################################
## Now customise build to comply with BUILD_FLAGS: ##
#####################################################

SOURCES = src/euler3d_cpu_double.cpp \
		  src/Base/common.cpp \
		  src/Base/config.cpp \
		  src/Base/io.cpp \
		  src/Base/io_enhanced.cpp \
		  src/Kernels/flux_kernels.cpp \
		  src/Kernels/kernels.cpp \
		  src/Kernels/mg.cpp \
		  src/Kernels/indirect_rw_kernel.cpp \
		  src/Kernels/validation.cpp \
		  src/Monitoring/timer.cpp \
		  src/Monitoring/papi_funcs.cpp \
		  src/Monitoring/loop_stats.cpp

OBJECTS     := $(patsubst src/%.cpp, $(OBJ_DIR)/%.o,     $(SOURCES))
OBJECTS_DBG := $(patsubst src/%.cpp, $(OBJ_DIR_DBG)/%.o, $(SOURCES))

INCLUDES := -Isrc -Isrc/Base -Isrc/Kernels -Isrc/Monitoring

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
debug: generic

$(OBJ_DIR)/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CPP) $(CFLAGS) $(OPTIMISATION) -c -o $@ $< $(BUILD_FLAGS) $(INCLUDES)

generic: $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(CPP) $(CFLAGS) $(OPTIMISATION) $^ -o $(BIN_NAME) $(BUILD_FLAGS) $(INCLUDES) $(LIBS)

clean:
	rm -rf $(BIN_DIR)
	rm -rf obj/$(shell hostname)
