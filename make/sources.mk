#===============
# Source files
#===============

STATIC_LIB := lib/libgemmul8.a
SHARED_LIB := lib/libgemmul8.so

ifeq ($(filter clean,$(MAKECMDGOALS)),)

HEADER_DIR := include src
HEADER := $(sort $(shell \
    for d in $(HEADER_DIR); do \
        [ -d $$d ] && find $$d -name '*.hpp'; \
    done))

CU_SRC_DIR := src
CU_SRC := $(sort $(shell find $(CU_SRC_DIR) -type f -path '*/cu/*.cu'))

HOOK_DIR := hook
HOOK_SRC := $(sort $(wildcard $(HOOK_DIR)/destroy.cu $(addprefix $(HOOK_DIR)/,$(addsuffix .cu,$(BUILD_OPS)))))

OBJ_DIR := build/$(BACKEND)/$(BUILD_VARIANT)
obj_from_src = $(OBJ_DIR)/$(patsubst %.cu,%.o,$(subst /cu/,/,$(1)))

CU_OBJ := $(foreach f,$(CU_SRC),$(call obj_from_src,$(f)))
HOOK_OBJ := $(foreach f,$(HOOK_SRC),$(call obj_from_src,$(f)))

endif
