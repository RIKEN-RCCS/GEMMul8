#===============
# GEMMul8 build
#===============

.DEFAULT_GOAL := all

include make/config.mk
include make/backend.mk
include make/functions.mk
include make/selection.mk
include make/sources.mk
include make/rules.mk
include make/inst_common.mk

export TMPDIR := $(TEMPDIR)

ifeq ($(filter clean,$(MAKECMDGOALS)),)
INST_MK := $(sort $(shell find src -name instantiations.mk 2>/dev/null))
include $(INST_MK)
ALL_OBJ := $(CU_OBJ) $(INST_OBJ)
else
ALL_OBJ :=
endif

LINK_FLAGS :=

ifeq ($(shell uname -m),aarch64)
LINK_FLAGS += -Xlinker --stub-group-size=33554432
endif

.PHONY: all clean info compile_objects_banner FORCE build-plan

all: info $(STATIC_LIB) $(SHARED_LIB)

$(STATIC_LIB): $(ALL_OBJ) $(ACTIVE_CONFIG) $(BUILD_MAKEFILES)
	@mkdir -p lib
	@echo "Creating static library"
	@echo "AR  $@"
	@rm -f $@.tmp
	@ar rcs $@.tmp $(ALL_OBJ)
	@mv -f $@.tmp $@

$(SHARED_LIB): $(ALL_OBJ) $(HOOK_OBJ) $(ACTIVE_CONFIG) $(BUILD_MAKEFILES)
	@mkdir -p lib
	@echo "Creating shared library"
	@echo "LD  $@"
	@$(COMPILER) $(ARCH) -shared -o $@.tmp $(ALL_OBJ) $(HOOK_OBJ) $(LIBS) $(LINK_FLAGS)
	@mv -f $@.tmp $@

$(ACTIVE_CONFIG): FORCE
	@mkdir -p $(@D)
	@printf '%s\n' '$(BUILD_SIGNATURE)' > $@.tmp
	@cmp -s $@.tmp $@ && rm -f $@.tmp || mv -f $@.tmp $@

FORCE:

build-plan:
	@echo 'OPS          : $(BUILD_OPS)'
	@echo 'OZ2_BACKENDS : $(BUILD_OZ2_BACKENDS)'
	@echo 'OBJECTS      : $(words $(ALL_OBJ)) + $(words $(HOOK_OBJ)) hooks'
	@$(foreach obj,$(ALL_OBJ) $(HOOK_OBJ),echo '$(obj)';)

clean:
	rm -rf build lib $(COMPILE_INFO)
