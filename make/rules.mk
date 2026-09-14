#===============
# User-visible build info
#===============

COMPILE_INFO := compile_info

HOST_CXX ?= $(CXX)
ifeq ($(strip $(HOST_CXX)),)
HOST_CXX := c++
endif

info:
	@{ \
	    echo ""; \
	    echo "BACKEND      : $(BACKEND)"; \
	    echo "OPS          : $(BUILD_OPS)"; \
	    echo "OZ2_BACKENDS : $(BUILD_OZ2_BACKENDS)"; \
	    if [ "$(BACKEND)" = "cuda" ]; then \
	        echo "CUDA_PATH    : $(CUDA_PATH)"; \
	    fi; \
	    if [ "$(BACKEND)" = "hip" ]; then \
	        echo "HIP_PATH     : $(HIP_PATH)"; \
	    fi; \
	    echo "COMPILER     : $(COMPILER)"; \
		host_cxx="$(HOST_CXX)"; \
	    host_cxx_path="$$(command -v "$$host_cxx" 2>/dev/null || true)"; \
	    [ -n "$$host_cxx_path" ] || host_cxx_path="$$host_cxx"; \
	    host_cxx_version="$$( \
	        "$$host_cxx" --version 2>/dev/null | head -n 1 \
	    )"; \
	    [ -n "$$host_cxx_version" ] || host_cxx_version="unknown"; \
	    echo "HOST_CXX     : $$host_cxx_path"; \
	    echo "HOST_CXX_VER : $$host_cxx_version"; \
	    echo "KERNEL_BUILD : $$kernel_build"; \
	    echo "GPU_ARCH     : $(GPU_ARCH)"; \
		cpu_model="$$( \
	        if command -v lscpu >/dev/null 2>&1; then \
	            lscpu | awk -F: '/^Model name:/ {gsub(/^[[:space:]]+/, "", $$2); print $$2; exit}'; \
	        elif [ -r /proc/cpuinfo ]; then \
	            awk -F: '/^model name[[:space:]]*:/ {gsub(/^[[:space:]]+/, "", $$2); print $$2; exit}' /proc/cpuinfo; \
	        elif command -v sysctl >/dev/null 2>&1; then \
	            sysctl -n machdep.cpu.brand_string 2>/dev/null; \
	        fi); \
	    "; \
	    cpu_threads="$$( \
	        if command -v nproc >/dev/null 2>&1; then \
	            nproc; \
	        elif command -v getconf >/dev/null 2>&1; then \
	            getconf _NPROCESSORS_ONLN 2>/dev/null; \
	        elif command -v sysctl >/dev/null 2>&1; then \
	            sysctl -n hw.logicalcpu 2>/dev/null; \
	        fi); \
	    "; \
	    cpu_cores="$$( \
	        if command -v lscpu >/dev/null 2>&1; then \
	            lscpu -p=CORE,SOCKET 2>/dev/null | awk -F, '!/^#/ {print $$1 "," $$2}' | sort -u | wc -l | tr -d ' '; \
	        elif command -v sysctl >/dev/null 2>&1; then \
	            sysctl -n hw.physicalcpu 2>/dev/null; \
	        fi); \
	    "; \
	    [ -n "$$cpu_model" ] || cpu_model="unknown"; \
	    [ -n "$$cpu_threads" ] || cpu_threads="unknown"; \
	    [ -n "$$cpu_cores" ] || cpu_cores="unknown"; \
	    echo "CPU          : $$cpu_model"; \
	    echo "CPU_CORES    : $$cpu_cores"; \
	    echo "CPU_THREADS  : $$cpu_threads"; \
	    if [ "$(BACKEND)" = "cuda" ]; then \
	        if command -v nvidia-smi >/dev/null 2>&1; then \
	            smi_info="$$(nvidia-smi --version 2>/dev/null)" && { \
	                echo ""; \
	                echo "$$smi_info"; \
	            }; \
	        fi; \
	    fi; \
	    if [ "$(BACKEND)" = "hip" ]; then \
	        if command -v amd-smi >/dev/null 2>&1; then \
	            smi_info="$$(amd-smi version 2>/dev/null)" && { \
	                echo ""; \
	                echo "$$smi_info"; \
	            }; \
	        fi; \
	    fi; \
	    echo ""; \
	} | tee $(COMPILE_INFO)

compile_objects_banner: info
	@mkdir -p $(TEMPDIR)
	@echo "Compiling objects"


#===============
# Normal .cu compile rule
#===============

ifeq ($(filter clean,$(MAKECMDGOALS)),)

define COMPILE_CU
$(call obj_from_src,$(1)): $(1) $(HEADER) $(BUILD_MAKEFILES) | compile_objects_banner
	@mkdir -p $$(dir $$@)
	@echo "$$(COMPILER)  -c $$<  -o $$@"
	@$$(COMPILER) $$(FLAGS_PIC) $$(ARCH) -c $$< -o $$@
endef

$(foreach f,$(CU_SRC) $(HOOK_SRC),$(eval $(call COMPILE_CU,$(f))))

endif
