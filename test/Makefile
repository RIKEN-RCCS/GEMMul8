#===============
# options
#===============
# Select GPU backend: 
#   cuda, hip, or auto (auto-detect using nvidia-smi/rocminfo)
BACKEND ?= auto

# Operation families (each includes its ordinary and Lt API).
# Example: OPS="gemm trsm". Default: every operation.
OPS ?= all

# Ozaki-II arithmetic backends, independent of BACKEND=cuda|hip.
# Example: OZ2_BACKENDS=INT8 or OZ2_BACKENDS=FP8.
OZ2_BACKENDS ?= INT8 FP8

# Path to CUDA toolkit installation
CUDA_PATH ?= /usr/local/cuda

# Path to HIP (ROCm) toolkit installation
HIP_PATH ?= /opt/rocm

# Target GPU architecture
# e.g., 
#   GPU_ARCH=90     ... build for H100
#   GPU_ARCH=gfx942 ... build for MI300X
#   GPU_ARCH=auto   ... auto-detect using nvidia-smi/amd-smi
# see https://developer.nvidia.com/cuda/gpus for NVIDIA GPUs
# see https://rocm.docs.amd.com/en/latest/reference/gpu-arch-specs.html for AMD GPUs
GPU_ARCH ?= auto


ALL_OPS := gemm symm syrk syr2k syrkx hemm herk her2k herkx trmm trsm trtrmm
BUILD_OPS := $(sort $(strip $(OPS)))
ifneq ($(filter all,$(BUILD_OPS)),)
BUILD_OPS := $(sort $(ALL_OPS))
endif
ifneq ($(filter-out all $(ALL_OPS),$(OPS)),)
$(error Unknown OPS: $(filter-out all $(ALL_OPS),$(OPS)))
endif
ifeq ($(BUILD_OPS),)
$(error OPS must contain at least one operation or all)
endif

BUILD_OZ2_BACKENDS := $(sort $(strip $(OZ2_BACKENDS)))
ifneq ($(filter-out INT8 FP8,$(BUILD_OZ2_BACKENDS)),)
$(error OZ2_BACKENDS must contain INT8 and/or FP8)
endif
ifeq ($(BUILD_OZ2_BACKENDS),)
$(error OZ2_BACKENDS must contain INT8 and/or FP8)
endif

BUILD_DEFINES := $(foreach b,INT8 FP8,-DGEMMUL8_BUILD_$(b)=$(if $(filter $(b),$(BUILD_OZ2_BACKENDS)),1,0))
BUILD_DEFINES += $(foreach op,$(ALL_OPS),-DGEMMUL8_BUILD_OP_$(op)=$(if $(filter $(op),$(BUILD_OPS)),1,0))

#===============
# Auto-detect backend (CUDA or HIP)
#===============

ifeq ($(BACKEND),auto)
ifneq ($(shell command -v nvidia-smi 2>/dev/null),)
override BACKEND := cuda
else ifneq ($(shell command -v rocminfo 2>/dev/null),)
override BACKEND := hip
else ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),run)
$(error Neither NVIDIA (CUDA) nor AMD (ROCm) GPU environment detected!)
endif
endif
endif


#===============
# CUDA setup
#===============

ifeq ($(BACKEND),cuda)

ifeq ($(GPU_ARCH),auto)
override GPU_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')
endif

export PATH := $(CUDA_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(CUDA_PATH)/lib64:$(LD_LIBRARY_PATH)

COMPILER := nvcc
FLAGS := -std=c++20 -O3 -diag-suppress 177 -I../include
LIBS := -lcublas -lcublasLt -lcudart -lcuda -lnvidia-ml -ldl ../lib/libgemmul8.a
ARCH := -gencode arch=compute_$(GPU_ARCH),code=sm_$(GPU_ARCH)

endif


#===============
# HIP setup
#===============

ifeq ($(BACKEND),hip)

ifeq ($(GPU_ARCH),auto)
override GPU_ARCH := $(shell amd-smi static --asic --csv | grep -oE 'gfx[0-9]+[a-z]*' | head -n 1)
endif

export PATH := $(HIP_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(HIP_PATH)/lib:$(LD_LIBRARY_PATH)

COMPILER := hipcc
FLAGS := -std=c++20 -O3 -I../include
FLAGS += -ffp-contract=off
FLAGS += -Wno-unused-result -Wno-unused-command-line-argument -Wno-unused-value
FLAGS += -DOCML_BASIC_ROUNDED_OPERATIONS
LIBS := -lamd_smi -lamdhip64 -lhipblas -lhipblaslt -ldl -Wl,../lib/libgemmul8.a
ARCH := --offload-arch=$(GPU_ARCH)

endif

FLAGS += $(BUILD_DEFINES)

#===============
# Files
#===============
TARGET := test

TIME_SRC := $(addsuffix .cu,$(addprefix time/time_,$(BUILD_OPS)))
ACCURACY_SRC := $(addsuffix .cu,$(addprefix accuracy/accuracy_,$(BUILD_OPS)))

CU_SRC := test.cu $(TIME_SRC) $(ACCURACY_SRC)
CU_OBJ := $(CU_SRC:.cu=.o)

TEMPDIR ?= tmp
export TMPDIR := $(abspath $(TEMPDIR))
COMPILE_INFO := compile_info
ACTIVE_CONFIG := $(TEMPDIR)/active-config

HOST_CXX ?= $(CXX)
ifeq ($(strip $(HOST_CXX)),)
HOST_CXX := c++
endif

#===============
# compile
#===============

.PHONY: all clean info run FORCE

all: $(TARGET)

$(TARGET): $(CU_OBJ)
	$(COMPILER) $(FLAGS) $(ARCH) $(CU_OBJ) -o $@ $(LIBS)

FORCE:

$(ACTIVE_CONFIG): FORCE | $(TEMPDIR)
	@{ \
	    echo "BACKEND=$(BACKEND)"; \
	    echo "GPU_ARCH=$(GPU_ARCH)"; \
	    echo "OPS=$(BUILD_OPS)"; \
	    echo "OZ2_BACKENDS=$(BUILD_OZ2_BACKENDS)"; \
	} > $@.tmp
	@if [ ! -f $@ ] || ! cmp -s $@.tmp $@; then \
	    mv $@.tmp $@; \
	else \
	    rm -f $@.tmp; \
	fi

%.o: %.cu Makefile $(ACTIVE_CONFIG) | $(TEMPDIR) info
	@mkdir -p $(dir $@)
	$(COMPILER) $(FLAGS) $(ARCH) -c $< -o $@

$(TEMPDIR):
	mkdir -p $(TEMPDIR)

info: | $(TEMPDIR)
	@{ \
	    echo ""; \
	    echo "BACKEND      : $(BACKEND)"; \
	    echo "OPS          : $(BUILD_OPS)"; \
	    echo "OZ2_BACKENDS : $(BUILD_OZ2_BACKENDS)"; \
	    if [ "$(BACKEND)" = "cuda" ]; then \
	        echo "CUDA_PATH    : $(CUDA_PATH)"; \
	    fi; \
	    if [ "$(BACKEND)" = "hip" ]; then \
	        echo "HIP_PATH     : $(HIP_PATH)"; \
	    fi; \
	    echo "COMPILER     : $(COMPILER)"; \
		host_cxx="$(HOST_CXX)"; \
	    host_cxx_path="$$(command -v "$$host_cxx" 2>/dev/null || true)"; \
	    [ -n "$$host_cxx_path" ] || host_cxx_path="$$host_cxx"; \
	    host_cxx_version="$$( \
	        "$$host_cxx" --version 2>/dev/null | head -n 1 \
	    )"; \
	    [ -n "$$host_cxx_version" ] || host_cxx_version="unknown"; \
	    echo "HOST_CXX     : $$host_cxx_path"; \
	    echo "HOST_CXX_VER : $$host_cxx_version"; \
	    echo "GPU_ARCH     : $(GPU_ARCH)"; \
	    cpu_model="$$( \
	        if command -v lscpu >/dev/null 2>&1; then \
	            lscpu | awk -F: '/^Model name:/ {gsub(/^[[:space:]]+/, "", $$2); print $$2; exit}'; \
	        elif [ -r /proc/cpuinfo ]; then \
	            awk -F: '/^model name[[:space:]]*:/ {gsub(/^[[:space:]]+/, "", $$2); print $$2; exit}' /proc/cpuinfo; \
	        elif command -v sysctl >/dev/null 2>&1; then \
	            sysctl -n machdep.cpu.brand_string 2>/dev/null; \
	        fi); \
	    "; \
	    cpu_threads="$$( \
	        if command -v nproc >/dev/null 2>&1; then \
	            nproc; \
	        elif command -v getconf >/dev/null 2>&1; then \
	            getconf _NPROCESSORS_ONLN 2>/dev/null; \
	        elif command -v sysctl >/dev/null 2>&1; then \
	            sysctl -n hw.logicalcpu 2>/dev/null; \
	        fi); \
	    "; \
	    cpu_cores="$$( \
	        if command -v lscpu >/dev/null 2>&1; then \
	            lscpu -p=CORE,SOCKET 2>/dev/null | awk -F, '!/^#/ {print $$1 "," $$2}' | sort -u | wc -l | tr -d ' '; \
	        elif command -v sysctl >/dev/null 2>&1; then \
	            sysctl -n hw.physicalcpu 2>/dev/null; \
	        fi); \
	    "; \
	    [ -n "$$cpu_model" ] || cpu_model="unknown"; \
	    [ -n "$$cpu_threads" ] || cpu_threads="unknown"; \
	    [ -n "$$cpu_cores" ] || cpu_cores="unknown"; \
	    echo "CPU          : $$cpu_model"; \
	    echo "CPU_CORES    : $$cpu_cores"; \
	    echo "CPU_THREADS  : $$cpu_threads"; \
	    if [ "$(BACKEND)" = "cuda" ]; then \
	        if command -v nvidia-smi >/dev/null 2>&1; then \
	            smi_info="$$(nvidia-smi --version 2>/dev/null)" && { \
	                echo ""; \
	                echo "$$smi_info"; \
	            }; \
	        fi; \
	    fi; \
	    if [ "$(BACKEND)" = "hip" ]; then \
	        if command -v amd-smi >/dev/null 2>&1; then \
	            smi_info="$$(amd-smi version 2>/dev/null)" && { \
	                echo ""; \
	                echo "$$smi_info"; \
	            }; \
	        fi; \
	    fi; \
	    echo ""; \
	} | tee $(COMPILE_INFO)

	

#===============
# run
#===============

run:
	./$(TARGET) $(MODE)


#===============
# Clean
#===============

clean:
	rm -f test.o
	rm -f time/*.o
	rm -f accuracy/*.o
	rm -f $(TARGET)
	rm -f $(COMPILE_INFO)
	rm -rf "$(TEMPDIR)"
