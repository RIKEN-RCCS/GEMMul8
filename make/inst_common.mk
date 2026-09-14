#===============
# Explicit-instantiation object rule generator
#===============

INST_OBJ :=

# $(1): object stem without build/$(BACKEND) and .o
# $(2): cu_recipe source file
# $(3): -D flags passed to nvcc/hipcc
ADD_INST_OBJ = $(if $(call inst_enabled,$(1),$(2),$(3)),$(call ADD_INST_OBJ_IMPL,$(1),$(2),$(3)))

define ADD_INST_OBJ_IMPL
INST_OBJ += $(OBJ_DIR)/$(strip $(1)).o

$(OBJ_DIR)/$(strip $(1)).o: $(strip $(2)) $(HEADER) $(BUILD_MAKEFILES) | compile_objects_banner
	@mkdir -p $$(dir $$@)
	@echo "$$(COMPILER)  -c $$<  -o $$@"
	@$$(COMPILER) $$(FLAGS_PIC) $$(ARCH) $(strip $(3)) -c $$< -o $$@
endef
