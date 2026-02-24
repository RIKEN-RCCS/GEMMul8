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
GPU_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')
endif

ifeq ($(GPU_MEM_MB),auto)
GPU_MEM_MB := $(shell nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -n 1)
endif

export PATH := $(CUDA_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(CUDA_PATH)/lib64:$(LD_LIBRARY_PATH)

COMPILER := nvcc
LIBS := -lcublas -lcublasLt -lcudart -lcuda -lnvidia-ml -ldl ../lib/libgemmul8.a
FLAGS := -std=c++20 -O3 -I../include
FLAGS += -DGPU_MEM_MB=$(GPU_MEM_MB)
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

ifeq ($(GPU_MEM_MB),auto)
GPU_MEM_MB := $(shell amd-smi static --vram --csv | awk -F, "NR==2 {print \$$4}")
endif

COMPILER := hipcc
LIBS := -lamd_smi -lamdhip64 -lhipblas -lhipblaslt -ldl -Wl,../lib/libgemmul8.a
FLAGS := -std=c++20 -O3 -I../include
FLAGS += -ffp-contract=off
FLAGS += -Wno-unused-result -Wno-unused-command-line-argument
FLAGS += -DOCML_BASIC_ROUNDED_OPERATIONS
FLAGS += -DGPU_MEM_MB=$(GPU_MEM_MB)
ARCH := --offload-arch=$(GPU_ARCH)

endif


#===============
# Compile
#===============

TARGET := test

all: INFO VERSION $(TARGET)

INFO:
	$(info BACKEND      : $(BACKEND))
ifeq ($(BACKEND),cuda)
	$(info CUDA_PATH    : $(CUDA_PATH))
endif
ifeq ($(BACKEND),hip)
	$(info HIP_PATH     : $(HIP_PATH))
endif
	$(info GPU_ARCH     : $(GPU_ARCH))
	$(info COMPILER     : $(COMPILER))
	$(info GPU_MEM [MB] : $(GPU_MEM_MB))

$(TARGET): $(TARGET).cu
	$(COMPILER) $< $(FLAGS) $(ARCH) -o $@ $(LIBS)

VERSION:
	$(COMPILER) --version 2>&1 | tee $(COMPILER)_version

run:
	./$(TARGET) $(MODE)

clean:
	rm -f *.o
	rm -f $(TARGET)
