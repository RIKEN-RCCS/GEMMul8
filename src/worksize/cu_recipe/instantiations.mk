#===============
# Explicit instantiations: src/worksize
#===============

WORKSIZE_COMPLEX_CODES := real complex
WORKSIZE_BACKEND_CODES := i8 f8

WORKSIZE_DIR := src/worksize
WORKSIZE_SRC := $(WORKSIZE_DIR)/cu_recipe/worksize.cu
SAVE_MEM_SRC := $(WORKSIZE_DIR)/cu_recipe/memory_saving.cu

# memory_saving feature
$(eval $(call ADD_INST_OBJ,\
    $(WORKSIZE_DIR)/memory_saving,\
    $(SAVE_MEM_SRC),\
    ))

# $(1): real|complex, $(2): i8|f8
define ADD_WORKSIZE
$(eval $(call ADD_INST_OBJ,\
    $(WORKSIZE_DIR)/worksize_$(1)_$(2),\
    $(WORKSIZE_SRC),\
    -DGEMMUL8_INST_COMPLEX=$(call complex_cpp_01,$(1)) \
    -DGEMMUL8_INST_BACKEND=$(call backend_cpp,$(2))))
endef

$(foreach cx,$(WORKSIZE_COMPLEX_CODES),\
  $(foreach b,$(WORKSIZE_BACKEND_CODES),\
    $(call ADD_WORKSIZE,$(cx),$(b))))
