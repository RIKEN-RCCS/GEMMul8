#===============
# Explicit instantiations: src/oz2/scaling/fast
#===============

FAST_TYPE_CODES_REAL := s d
FAST_TYPE_CODES_COMPLEX := c z
FAST_TYPE_CODES_ALL := $(FAST_TYPE_CODES_REAL) $(FAST_TYPE_CODES_COMPLEX)
FAST_BACKEND_CODES := i8 f8
FAST_TRI_PATTERNS := full_nonunit upper_nonunit upper_unit lower_nonunit lower_unit
FAST_UPLO_CODES := upper lower

FAST_DIR := src/oz2/scaling/fast
FAST_SCALING_SRC := $(FAST_DIR)/cu_recipe/scaling.cu
FAST_SYMM_SRC    := $(FAST_DIR)/cu_recipe/scaling_symm.cu
FAST_HEMM_SRC    := $(FAST_DIR)/cu_recipe/scaling_hemm.cu

# $(1): s|d|c|z, $(2): i8|f8, $(3): triangular pattern
define ADD_FAST_SCALING
$(eval $(call ADD_INST_OBJ,\
    $(FAST_DIR)/scaling_$(1)_$(2)_$(3),\
    $(FAST_SCALING_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(call tri_uplo,$(3))) \
    -DGEMMUL8_INST_DIAG=$(call diag_cpp,$(call tri_diag,$(3)))))
endef

$(foreach t,$(FAST_TYPE_CODES_ALL),\
  $(foreach b,$(FAST_BACKEND_CODES),\
    $(foreach p,$(FAST_TRI_PATTERNS),\
      $(call ADD_FAST_SCALING,$(t),$(b),$(p)))))

# $(1): s|d|c|z, $(2): i8|f8, $(3): upper|lower
define ADD_FAST_SYMM
$(eval $(call ADD_INST_OBJ,\
    $(FAST_DIR)/scaling_symm_$(1)_$(2)_$(3),\
    $(FAST_SYMM_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach t,$(FAST_TYPE_CODES_ALL),\
  $(foreach b,$(FAST_BACKEND_CODES),\
    $(foreach u,$(FAST_UPLO_CODES),\
      $(call ADD_FAST_SYMM,$(t),$(b),$(u)))))

# $(1): c|z, $(2): i8|f8, $(3): upper|lower
define ADD_FAST_HEMM
$(eval $(call ADD_INST_OBJ,\
    $(FAST_DIR)/scaling_hemm_$(1)_$(2)_$(3),\
    $(FAST_HEMM_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach t,$(FAST_TYPE_CODES_COMPLEX),\
  $(foreach b,$(FAST_BACKEND_CODES),\
    $(foreach u,$(FAST_UPLO_CODES),\
      $(call ADD_FAST_HEMM,$(t),$(b),$(u)))))

# HERK-specific FP8 complex input plans.
ifneq ($(filter herk,$(BUILD_OPS)),)
define ADD_FAST_HERK
$(eval $(call ADD_INST_OBJ,\
    $(FAST_DIR)/scaling_herk_$(1)_f8,\
    $(FAST_SCALING_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=FP8 \
    -DGEMMUL8_INST_FILLMODE=CUBLAS_FILL_MODE_FULL \
    -DGEMMUL8_INST_DIAG=CUBLAS_DIAG_NON_UNIT \
    -DGEMMUL8_INST_CONJ=false -DGEMMUL8_INST_HERK=true))
endef
$(foreach t,c z,$(call ADD_FAST_HERK,$(t)))
endif
