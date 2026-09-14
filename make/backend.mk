#===============
# Auto-detect backend (CUDA or HIP)
#===============

ifeq ($(BACKEND),auto)

ifneq ($(shell command -v nvidia-smi 2>/dev/null),)
override BACKEND := cuda

else ifneq ($(shell command -v rocminfo 2>/dev/null),)
override BACKEND := hip

else ifeq ($(filter clean,$(MAKECMDGOALS)),)
$(error Neither NVIDIA CUDA nor AMD ROCm environment detected)

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
FLAGS := -std=c++20 -O3 -diag-suppress 177 -diag-suppress 186 -Iinclude -Isrc -DGPU_ARCH=$(GPU_ARCH)
LIBS := -lcublas -lcudart -lcublasLt -lcuda -lnvidia-ml -ldl
ARCH := -gencode arch=compute_$(GPU_ARCH),code=sm_$(GPU_ARCH)
FLAGS_PIC := $(FLAGS) -Xcompiler -fPIC

endif


#===============
# HIP setup
#===============

ifeq ($(BACKEND),hip)

ifeq ($(GPU_ARCH),auto)
override GPU_ARCH := $(shell \
    if command -v amd-smi >/dev/null 2>&1; then \
        amd-smi static --asic --csv 2>/dev/null | grep -oE 'gfx[0-9]+[a-z]*' | head -n 1; \
    else \
        rocminfo 2>/dev/null | awk '/Name:[[:space:]]*gfx/ {print $$2; exit}'; \
    fi)
endif

export PATH := $(HIP_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(HIP_PATH)/lib:$(LD_LIBRARY_PATH)

COMPILER := hipcc
FLAGS := -std=c++20 -O3 -Iinclude -Isrc
FLAGS += -ffp-contract=off
FLAGS += -Wno-unused-result -Wno-unused-command-line-argument -Wno-unused-value
FLAGS += -DOCML_BASIC_ROUNDED_OPERATIONS -DGPU_ARCH=$(GPU_ARCH)
LIBS := -lamd_smi -lamdhip64 -lhipblas -lhipblaslt -ldl
ARCH := --offload-arch=$(GPU_ARCH)
FLAGS_PIC := $(FLAGS) -fPIC

endif
