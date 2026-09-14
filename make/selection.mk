# Operation / arithmetic selection and the shared-kernel dependency closure.
# Ordinary APIs are INT8-only; Lt APIs are emitted for every selected backend.

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

empty :=
space := $(empty) $(empty)
OPS_ID := $(if $(filter-out $(BUILD_OPS),$(ALL_OPS)),$(subst $(space),+,$(BUILD_OPS)),all)
OZ2_ID := $(subst $(space),+,$(BUILD_OZ2_BACKENDS))
BUILD_VARIANT := $(GPU_ARCH)/$(OPS_ID)/$(OZ2_ID)
BUILD_SIGNATURE := $(BACKEND)/$(BUILD_VARIANT)
ACTIVE_CONFIG := build/active-config
BUILD_MAKEFILES := Makefile $(wildcard make/*.mk) $(wildcard src/*/cu_recipe/instantiations.mk) \
                  $(wildcard src/oz2/*/cu_recipe/instantiations.mk) \
                  $(wildcard src/oz2/scaling/*/cu_recipe/instantiations.mk)

BUILD_DEFINES := $(foreach b,INT8 FP8,-DGEMMUL8_BUILD_$(b)=$(if $(filter $(b),$(BUILD_OZ2_BACKENDS)),1,0))
BUILD_DEFINES += $(foreach op,$(ALL_OPS),-DGEMMUL8_BUILD_OP_$(op)=$(if $(filter $(op),$(BUILD_OPS)),1,0))
FLAGS_PIC += $(BUILD_DEFINES)

# Every operation can use full-matrix scaling (including blocked off-diagonal
# products and the GEMM updates inside TRSM). TRSM's triangular solves are
# native BLAS calls, so they do not require triangular Ozaki scaling.
CORE_SOURCES := \
    src/oz2/scaling/fast/cu_recipe/scaling.cu \
    src/oz2/scaling/accu/cu_recipe/extract.cu \
    src/oz2/scaling/accu/cu_recipe/scaling.cu \
    src/oz2/scaling/general/cu_recipe/scaling_rowwise.cu \
    src/oz2/mod/cu_recipe/mod_hi2mid.cu \
    src/oz2/mod/cu_recipe/mod_reduce_matprod.cu \
    src/oz2/undo_scaling/cu_recipe/undo_scaling.cu

CORE_SOURCES += $(if $(filter symm,$(BUILD_OPS)),\
    src/oz2/scaling/fast/cu_recipe/scaling_symm.cu \
    src/oz2/scaling/accu/cu_recipe/extract_symm.cu \
    src/oz2/scaling/accu/cu_recipe/scaling_symm.cu \
    src/oz2/scaling/general/cu_recipe/scaling_symm.cu)
CORE_SOURCES += $(if $(filter hemm,$(BUILD_OPS)),\
    src/oz2/scaling/fast/cu_recipe/scaling_hemm.cu \
    src/oz2/scaling/accu/cu_recipe/extract_hemm.cu \
    src/oz2/scaling/accu/cu_recipe/scaling_hemm.cu \
    src/oz2/scaling/general/cu_recipe/scaling_hemm.cu)
CORE_SOURCES += $(if $(filter syrk,$(BUILD_OPS)),src/oz2/scaling/accu/cu_recipe/scaling_syrk.cu)
CORE_SOURCES += $(if $(filter herk,$(BUILD_OPS)),\
    src/oz2/scaling/accu/cu_recipe/scaling_herk.cu src/oz2/mod/cu_recipe/mod_hi2mid_aha.cu)
CORE_SOURCES += $(if $(filter syr2k,$(BUILD_OPS)),src/oz2/undo_scaling/cu_recipe/undo_scaling_syr2k.cu)
CORE_SOURCES += $(if $(filter her2k,$(BUILD_OPS)),src/oz2/undo_scaling/cu_recipe/undo_scaling_her2k.cu)

# Keep all workspace-query APIs for the selected arithmetic backends. They
# contain size arithmetic only, and are also used by hook workspace planning.
# Configuration APIs remain available regardless of OPS.
SELECTED_SOURCES := $(CORE_SOURCES) src/worksize/cu_recipe/%.cu \
    src/trsm/cu_recipe/block_size_trsm.cu \
    $(addsuffix /cu_recipe/%.cu,$(addprefix src/,$(BUILD_OPS)))

CORE_INPUT_UPLOS := CUBLAS_FILL_MODE_FULL \
    $(if $(filter trmm trtrmm,$(BUILD_OPS)),CUBLAS_FILL_MODE_UPPER CUBLAS_FILL_MODE_LOWER)
CORE_ACCU_OUTPUT_UPLOS := CUBLAS_FILL_MODE_FULL \
    $(if $(filter syrkx herkx trtrmm,$(BUILD_OPS)),CUBLAS_FILL_MODE_UPPER CUBLAS_FILL_MODE_LOWER)
CORE_OUTPUT_UPLOS := CUBLAS_FILL_MODE_FULL \
    $(if $(filter syrk herk syr2k her2k syrkx herkx trtrmm,$(BUILD_OPS)),CUBLAS_FILL_MODE_UPPER CUBLAS_FILL_MODE_LOWER)
CORE_TYPES := cuFloatComplex cuDoubleComplex \
    $(if $(filter-out hemm herk her2k herkx,$(BUILD_OPS)),float double)
CORE_UNDO_TRIPLES := float_float_float double_double_double \
    cuFloatComplex_cuFloatComplex_cuFloatComplex cuDoubleComplex_cuDoubleComplex_cuDoubleComplex \
    $(if $(filter herk,$(BUILD_OPS)),cuFloatComplex_float_float cuDoubleComplex_double_double) \
    $(if $(filter herkx,$(BUILD_OPS)),cuFloatComplex_cuFloatComplex_float cuDoubleComplex_cuDoubleComplex_double)

GENERIC_SCALING_SOURCES := src/oz2/scaling/fast/cu_recipe/scaling.cu \
    src/oz2/scaling/accu/cu_recipe/extract.cu src/oz2/scaling/accu/cu_recipe/scaling.cu \
    src/oz2/scaling/general/cu_recipe/scaling_rowwise.cu

# Inspect the recipe's explicit -D arguments, rather than object-name suffixes.
inst_value = $(patsubst -D$(1)=%,%,$(filter -D$(1)=%,$(2)))
inst_backend_ok = $(if $(call inst_value,GEMMUL8_INST_BACKEND,$(1)),\
    $(filter $(BUILD_OZ2_BACKENDS),$(call inst_value,GEMMUL8_INST_BACKEND,$(1))),yes)
inst_type_ok = $(if $(filter src/oz2/%,$(1)),\
    $(if $(call inst_value,GEMMUL8_INST_TYPE,$(2)),\
        $(filter $(CORE_TYPES),$(call inst_value,GEMMUL8_INST_TYPE,$(2))),yes),yes)
inst_input_ok = $(if $(filter $(GENERIC_SCALING_SOURCES),$(1)),\
    $(filter $(CORE_INPUT_UPLOS),$(call inst_value,GEMMUL8_INST_FILLMODE,$(2))),yes)
inst_output_ok = $(if $(filter src/oz2/scaling/accu/cu_recipe/scaling.cu,$(1)),\
    $(filter $(CORE_ACCU_OUTPUT_UPLOS),$(call inst_value,GEMMUL8_INST_FILLMODE_C,$(2))),\
    $(if $(filter src/oz2/mod/cu_recipe/mod_hi2mid.cu src/oz2/undo_scaling/cu_recipe/undo_scaling.cu,$(1)),\
        $(filter $(CORE_OUTPUT_UPLOS),$(call inst_value,GEMMUL8_INST_FILLMODE,$(2))),yes))
inst_undo_ok = $(if $(filter src/oz2/undo_scaling/cu_recipe/undo_scaling.cu,$(1)),\
    $(filter $(CORE_UNDO_TRIPLES),$(call inst_value,GEMMUL8_INST_TYPE,$(2))_$(call inst_value,GEMMUL8_INST_TYPE_ALPHA,$(2))_$(call inst_value,GEMMUL8_INST_TYPE_BETA,$(2))),yes)
inst_enabled = $(and $(filter $(SELECTED_SOURCES),$(strip $(2))),\
    $(strip $(call inst_backend_ok,$(3))),$(strip $(call inst_type_ok,$(2),$(3))),\
    $(strip $(call inst_input_ok,$(2),$(3))),$(strip $(call inst_output_ok,$(2),$(3))),\
    $(strip $(call inst_undo_ok,$(2),$(3))))
