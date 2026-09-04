#===============
# Explicit instantiations: src/oz2/mod
#===============

MOD_COMPLEX_CODES := real complex
MOD_BACKEND_CODES := i8 f8
MOD_UPLO_CODES := upper lower
MOD_UPLO_CODES_WITH_FULL := full upper lower

MOD_DIR := src/oz2/mod
MOD_HI2MID_SRC := $(MOD_DIR)/cu_recipe/mod_hi2mid.cu
MOD_HI2MID_AHA_SRC := $(MOD_DIR)/cu_recipe/mod_hi2mid_aha.cu
MOD_REDUCE_MATPROD_SRC := $(MOD_DIR)/cu_recipe/mod_reduce_matprod.cu

# $(1): real|complex, $(2): i8|f8, $(3): full|upper|lower
define ADD_MOD_HI2MID
$(eval $(call ADD_INST_OBJ,\
    $(MOD_DIR)/mod_hi2mid_$(1)_$(2)_$(3),\
    $(MOD_HI2MID_SRC),\
    -DGEMMUL8_INST_COMPLEX=$(call complex_cpp,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach cx,$(MOD_COMPLEX_CODES),\
  $(foreach b,$(MOD_BACKEND_CODES),\
    $(foreach u,$(MOD_UPLO_CODES_WITH_FULL),\
      $(call ADD_MOD_HI2MID,$(cx),$(b),$(u)))))

# $(1): i8|f8, $(2): upper|lower
define ADD_MOD_HI2MID_AHA
$(eval $(call ADD_INST_OBJ,\
    $(MOD_DIR)/mod_hi2mid_aha_$(1)_$(2),\
    $(MOD_HI2MID_AHA_SRC),\
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(1)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(2))))
endef

$(foreach b,$(MOD_BACKEND_CODES),\
  $(foreach u,$(MOD_UPLO_CODES),\
    $(call ADD_MOD_HI2MID_AHA,$(b),$(u))))

# $(1): i8|f8
define ADD_MOD_REDUCE_MATPROD
$(eval $(call ADD_INST_OBJ,\
    $(MOD_DIR)/mod_reduce_matprod_$(1),\
    $(MOD_REDUCE_MATPROD_SRC),\
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(1))))
endef

$(foreach b,$(MOD_BACKEND_CODES),\
  $(call ADD_MOD_REDUCE_MATPROD,$(b)))
