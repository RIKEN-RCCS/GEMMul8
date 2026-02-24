#===============
# options
#===============

CUDA_PATH ?= /usr/local/cuda
HIP_PATH ?= /opt/rocm
BACKEND ?= auto
GPU_ARCH ?= auto
GPU_MEM_MB ?= auto


#===============
# Auto-detect backend (CUDA or HIP)
#===============

ifeq ($(BACKEND),auto)
ifneq ($(shell command -v nvidia-smi 2>/dev/null),)
BACKEND := cuda
else ifneq ($(shell command -v rocminfo 2>/dev/null),)
BACKEND := hip
else
$(error Neither NVIDIA (CUDA) nor AMD (ROCm) GPU environment detected!)
endif
endif


#===============
# CUDA setup
#===============

ifeq ($(BACKEND),cuda)

ifeq ($(GPU_ARCH),auto)
GPU_ARCH = $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')
endif

export PATH := $(CUDA_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(CUDA_PATH)/lib64:$(LD_LIBRARY_PATH)

COMPILER := nvcc
FLAGS := -std=c++20 -O3 -diag-suppress 177
LIBS := -lcublas -lcudart -lcublasLt -lcuda -lnvidia-ml -ldl
ARCH := -gencode arch=compute_$(GPU_ARCH),code=sm_$(GPU_ARCH)

endif


#===============
# HIP setup
#===============

ifeq ($(BACKEND),hip)

ifeq ($(GPU_ARCH),auto)
GPU_ARCH := $(shell amd-smi static --asic --csv | grep -o 'gfx[0-9]\+' | head -n 1)
endif

export PATH := $(HIP_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(HIP_PATH)/lib:$(LD_LIBRARY_PATH)

COMPILER := hipcc
FLAGS := -std=c++20 -O3
FLAGS += -ffp-contract=off
FLAGS += -Wno-unused-result -Wno-unused-command-line-argument
FLAGS += -DOCML_BASIC_ROUNDED_OPERATIONS
LIBS := -lamd_smi -lamdhip64 -lhipblas -lhipblaslt -ldl
ARCH := --offload-arch=$(GPU_ARCH)

endif


#===============
# Common
#===============

FLAGS_PIC := $(FLAGS) -Xcompiler -fPIC

MAKE_FLAGS := \
CUDA_PATH=$(CUDA_PATH) \
HIP_PATH=$(HIP_PATH) \
BACKEND=$(BACKEND) \
GPU_ARCH=$(GPU_ARCH) \
GPU_MEM_MB=$(GPU_MEM_MB)


#===============
# Files
#===============

HEADERS := include/gemmul8.hpp
CU_SRC := src/gemmul8.cu
CU_OBJ := $(CU_SRC:.cu=.o)
HOOK_SRC := src/hook.cu
HOOK_OBJ := src/hook.o
STATIC_LIB := lib/libgemmul8.a
SHARED_LIB := lib/libgemmul8.so


#===============
# Targets
#===============

.PHONY: all clean testing

all: $(STATIC_LIB) $(SHARED_LIB) testing


#===============
# Compile
#===============

%.o: %.cu $(HEADERS)
	$(COMPILER) $(FLAGS_PIC) $(ARCH) -c $< -o $@


#===============
# Static library
#===============

$(STATIC_LIB): $(CU_OBJ)
	mkdir -p lib
	ar rcs $@ $^


#===============
# Shared library
#===============

$(SHARED_LIB): $(CU_OBJ) $(HOOK_OBJ)
	mkdir -p lib
	$(COMPILER) -shared -o $@ $^ $(LIBS)


#===============
# Testing
#===============

testing: $(STATIC_LIB)
	$(MAKE) -C testing $(MAKE_FLAGS)


#===============
# Clean
#===============

clean:
	$(MAKE) -C testing clean
	rm -f \
	    $(CU_OBJ) \
	    $(HOOK_OBJ) \
	    $(STATIC_LIB) \
	    $(SHARED_LIB)
