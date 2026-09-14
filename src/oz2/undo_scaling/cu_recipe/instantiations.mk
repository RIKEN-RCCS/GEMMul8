#===============
# Explicit instantiations: src/oz2/undo_scaling
#===============

UNDO_BACKEND_CODES := i8 f8
UNDO_UPLO_CODES := upper lower
UNDO_UPLO_CODES_WITH_FULL := full upper lower
UNDO_TRIPLES := s_s_s d_d_d c_c_c c_c_s c_s_s z_z_z z_z_d z_d_d
UNDO_SYR2K_TRIPLES := s_s_s d_d_d c_c_c z_z_z
UNDO_HER2K_TRIPLES := c_c_s z_z_d

UNDO_DIR := src/oz2/undo_scaling
UNDO_SRC       := $(UNDO_DIR)/cu_recipe/undo_scaling.cu
UNDO_SYR2K_SRC := $(UNDO_DIR)/cu_recipe/undo_scaling_syr2k.cu
UNDO_HER2K_SRC := $(UNDO_DIR)/cu_recipe/undo_scaling_her2k.cu

# $(1): triplet T_TALPHA_TBETA, $(2): i8|f8, $(3): full|upper|lower
define ADD_UNDO
$(eval $(call ADD_INST_OBJ,\
    $(UNDO_DIR)/undo_scaling_$(call trip_T,$(1))_$(2)_$(3)_$(call trip_A,$(1))_$(call trip_B,$(1)),\
    $(UNDO_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(call trip_T,$(1))) \
    -DGEMMUL8_INST_TYPE_ALPHA=$(call type_cpp,$(call trip_A,$(1))) \
    -DGEMMUL8_INST_TYPE_BETA=$(call type_cpp,$(call trip_B,$(1))) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3)) \
    -DGEMMUL8_INST_HERMITIAN=$(if $(and $(filter c_c_s c_s_s z_z_d z_d_d,$(1)),$(filter upper lower,$(3))),1,0)))
endef

$(foreach trip,$(UNDO_TRIPLES),\
  $(foreach b,$(UNDO_BACKEND_CODES),\
    $(foreach u,$(UNDO_UPLO_CODES_WITH_FULL),\
      $(call ADD_UNDO,$(trip),$(b),$(u)))))

# $(1): triplet T_TALPHA_TBETA, $(2): i8|f8, $(3): upper|lower
define ADD_UNDO_SYR2K
$(eval $(call ADD_INST_OBJ,\
    $(UNDO_DIR)/undo_scaling_syr2k_$(call trip_T,$(1))_$(2)_$(3)_$(call trip_A,$(1))_$(call trip_B,$(1)),\
    $(UNDO_SYR2K_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(call trip_T,$(1))) \
    -DGEMMUL8_INST_TYPE_ALPHA=$(call type_cpp,$(call trip_A,$(1))) \
    -DGEMMUL8_INST_TYPE_BETA=$(call type_cpp,$(call trip_B,$(1))) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach trip,$(UNDO_SYR2K_TRIPLES),\
  $(foreach b,$(UNDO_BACKEND_CODES),\
    $(foreach u,$(UNDO_UPLO_CODES),\
      $(call ADD_UNDO_SYR2K,$(trip),$(b),$(u)))))

# $(1): triplet T_TALPHA_TBETA, $(2): i8|f8, $(3): upper|lower
define ADD_UNDO_HER2K
$(eval $(call ADD_INST_OBJ,\
    $(UNDO_DIR)/undo_scaling_her2k_$(call trip_T,$(1))_$(2)_$(3)_$(call trip_A,$(1))_$(call trip_B,$(1)),\
    $(UNDO_HER2K_SRC),\
    -DGEMMUL8_INST_TYPE=$(call type_cpp,$(call trip_T,$(1))) \
    -DGEMMUL8_INST_TYPE_ALPHA=$(call type_cpp,$(call trip_A,$(1))) \
    -DGEMMUL8_INST_TYPE_BETA=$(call type_cpp,$(call trip_B,$(1))) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2)) \
    -DGEMMUL8_INST_FILLMODE=$(call uplo_cpp,$(3))))
endef

$(foreach trip,$(UNDO_HER2K_TRIPLES),\
  $(foreach b,$(UNDO_BACKEND_CODES),\
    $(foreach u,$(UNDO_UPLO_CODES),\
      $(call ADD_UNDO_HER2K,$(trip),$(b),$(u)))))
