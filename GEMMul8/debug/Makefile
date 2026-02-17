
GPU_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')

CUDA_PATH ?= /usr/local/cuda
export PATH := $(CUDA_PATH)/bin:$(PATH)
export LD_LIBRARY_PATH := $(CUDA_PATH)/lib64:$(LD_LIBRARY_PATH)

COMPILER := nvcc
LIBS := -lcublas -lcudart -lcublasLt -lcuda -lnvidia-ml -ldl ../lib/libgemmul8.a
FLAGS := -std=c++20 -O3 -I../include
ARCH := -gencode arch=compute_$(GPU_ARCH),code=sm_$(GPU_ARCH)

export GEMMUL8_BACKEND=INT8
export GEMMUL8_NUM_MOD_D=15
export GEMMUL8_NUM_MOD_S=7
export GEMMUL8_NUM_MOD_Z=15
export GEMMUL8_NUM_MOD_C=7
export GEMMUL8_FASTMODE_D=0
export GEMMUL8_FASTMODE_S=1
export GEMMUL8_FASTMODE_Z=0
export GEMMUL8_FASTMODE_C=1
export GEMMUL8_MAXWS_BACKEND=INT8
export GEMMUL8_MAX_M=4096
export GEMMUL8_MAX_N=4096
export GEMMUL8_MAX_K=4096
export GEMMUL8_MAX_NUM_MOD=18
export GEMMUL8_SKIP_SCALE_A=1
export GEMMUL8_SKIP_SCALE_B=1

CFLAGS ?=
# CFLAGS += -DSGEMM
# CFLAGS += -DDGEMM
# CFLAGS += -DCGEMM
# CFLAGS += -DZGEMM
# CFLAGS += -DLT

LD:=../lib/libgemmul8.so
TARGET1 := test
TARGET2 := test_hijack

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(TARGET1).cu
	$(COMPILER) $< $(FLAGS) $(CFLAGS) $(ARCH) -o $@ $(LIBS)

$(TARGET2): $(TARGET2).cu
	$(COMPILER) $< $(FLAGS) $(ARCH) -o $@ $(LIBS)

run-test:
	compute-sanitizer --tool memcheck ./$(TARGET1)

run-test-hijack:
	LD_LIBRARY_PATH=$(CUDA_PATH)/lib64:$(LD_LIBRARY_PATH) LD_PRELOAD=$(LD):$(LD_PRELOAD) compute-sanitizer --tool memcheck ./$(TARGET2)

clean:
	rm -f *.o
	rm -f $(TARGET1) $(TARGET2)
