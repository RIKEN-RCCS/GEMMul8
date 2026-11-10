#===============
# Explicit instantiations: src/oz2/scaling/general
#===============

GENERAL_TYPE_CODES_REAL := s d
GENERAL_TYPE_CODES_COMPLEX := c z
GENERAL_TYPE_CODES_ALL := $(GENERAL_TYPE_CODES_REAL) $(GENERAL_TYPE_CODES_COMPLEX)
GENERAL_BACKEND_CODES := i8 f8
GENERAL_TRI_PATTERNS := full_nonunit upper_nonunit upper_unit lower_nonunit lower_unit
GENERAL_UPLO_CODES := upper lower
GENERAL_CONJ_CODES := conj0 conj1

GENERAL_DIR := src/oz2/scaling/general
GENERAL_ROWWISE_SRC := $(GENERAL_DIR)/cu_recipe/scaling_rowwise.cu
GENERAL_SYMM_SRC    := $(GENERAL_DIR)/cu_recipe/scaling_symm.cu
GENERAL_HEMM_SRC    := $(GENERAL_DIR)/cu_recipe/scaling_hemm.cu

# $(1): s|d|c|z, $(2): i8|f8, $(3): triangular pattern, $(4): conj0|conj1
define ADD_GENERAL_ROWWISE
$(eval $(call ADD_INST_OBJ,\
    $(GENERAL_DIR)/scaling_rowwise_$(1)_$(2)_$(3)_$(4),\
    $(GENERAL_ROWWISE_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(call tri_uplo,$(3))) \
    -DGEMMUL8_INST_DIAG=$(call diag_cpp,$(call tri_diag,$(3))) \
    -DGEMMUL8_INST_CONJ=$(call conj_cpp,$(4))))
endef

$(foreach t,$(GENERAL_TYPE_CODES_REAL),\
  $(foreach b,$(GENERAL_BACKEND_CODES),\
    $(foreach p,$(GENERAL_TRI_PATTERNS),\
      $(call ADD_GENERAL_ROWWISE,$(t),$(b),$(p),conj0))))

$(foreach t,$(GENERAL_TYPE_CODES_COMPLEX),\
  $(foreach b,$(GENERAL_BACKEND_CODES),\
    $(foreach p,$(GENERAL_TRI_PATTERNS),\
      $(foreach cj,$(GENERAL_CONJ_CODES),\
        $(call ADD_GENERAL_ROWWISE,$(t),$(b),$(p),$(cj))))))

# $(1): s|d|c|z, $(2): i8|f8, $(3): upper|lower
define ADD_GENERAL_SYMM
$(eval $(call ADD_INST_OBJ,\
    $(GENERAL_DIR)/scaling_symm_$(1)_$(2)_$(3),\
    $(GENERAL_SYMM_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach t,$(GENERAL_TYPE_CODES_ALL),\
  $(foreach b,$(GENERAL_BACKEND_CODES),\
    $(foreach u,$(GENERAL_UPLO_CODES),\
      $(call ADD_GENERAL_SYMM,$(t),$(b),$(u)))))

# $(1): c|z, $(2): i8|f8, $(3): upper|lower
define ADD_GENERAL_HEMM
$(eval $(call ADD_INST_OBJ,\
    $(GENERAL_DIR)/scaling_hemm_$(1)_$(2)_$(3),\
    $(GENERAL_HEMM_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach t,$(GENERAL_TYPE_CODES_COMPLEX),\
  $(foreach b,$(GENERAL_BACKEND_CODES),\
    $(foreach u,$(GENERAL_UPLO_CODES),\
      $(call ADD_GENERAL_HEMM,$(t),$(b),$(u)))))

# HERK-specific FP8 complex input plans.
ifneq ($(filter herk,$(BUILD_OPS)),)
define ADD_GENERAL_HERK
$(eval $(call ADD_INST_OBJ,\
    $(GENERAL_DIR)/scaling_herk_$(1)_f8_$(2),\
    $(GENERAL_ROWWISE_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=FP8 \
    -DGEMMUL8_INST_FILLMODE=CUBLAS_FILL_MODE_FULL \
    -DGEMMUL8_INST_DIAG=CUBLAS_DIAG_NON_UNIT \
    -DGEMMUL8_INST_CONJ=$(call conj_cpp,$(2)) -DGEMMUL8_INST_HERK=true))
endef
$(foreach t,c z,$(foreach cj,$(GENERAL_CONJ_CODES),$(call ADD_GENERAL_HERK,$(t),$(cj))))
endif
