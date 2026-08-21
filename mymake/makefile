CONFIG ?= m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)
include m_files

.DEFAULT_GOAL := all
.NOTPARALLEL:

obj_path := _OBJ
lib_path := LIB
library_name := libads.a
library_target := $(lib_path)/$(library_name)

MODULE_OUTPUT ?= -J.
FFLAGS += $(MODULE_OUTPUT) -I. -I./$(obj_path)
ARFLAGS = rcs

EXEC ?= l2
SOURCE_ALL ?= $(SOURCES) $(L2)
VPATH = $(dir $(SOURCE_ALL))

dot_o := .o
objF90 := $(notdir $(filter %$(dot_o),$(SOURCE_ALL:.F90=$(dot_o))))
OBJECT_ALL := $(addprefix $(obj_path)/,$(objF90))
core_objF90 := $(notdir $(filter %$(dot_o),$(SOURCES:.F90=$(dot_o))))
CORE_OBJECTS := $(addprefix $(obj_path)/,$(core_objF90))

PROBLEMS := l2 heat eriksson pure_diffusion_igrm oil igrm_l2 igrm_heat \
	igrm_eirksson

.PHONY: all build library problems list-problems clean clean-all \
	clean-problems clean-library clean-library-build help $(PROBLEMS)

# Backward-compatible low-level build: callers may still supply SOURCE_ALL and
# EXEC explicitly. The top-level Makefile uses the named clean-build targets.
all: build

build: EXEC/$(EXEC)

library:
	+$(MAKE) --no-print-directory -j1 clean-library-build \
		CONFIG="$(CONFIG_PATH)"
	+$(MAKE) --no-print-directory -j1 $(library_target) \
		CONFIG="$(CONFIG_PATH)"

EXEC/$(EXEC): $(OBJECT_ALL) | EXEC
	$(FF) $(FFLAGS) -o $@ $(OBJECT_ALL) $(USER_LIB)
	@printf 'Built %s\n' '$@'

$(library_target): $(CORE_OBJECTS) | $(lib_path)
	$(AR) $(ARFLAGS) $@ $(CORE_OBJECTS)
	@printf 'Built %s\n' '$@'

$(obj_path)/%$(dot_o): %.F90 $(MAKEFILE_LIST) | $(obj_path)
	$(FF) $(FFLAGS) -o $@ -c $<

EXEC $(obj_path) $(lib_path):
	mkdir -p $@

l2: PROBLEM_SOURCES := $(L2)
heat: PROBLEM_SOURCES := $(HEAT)
eriksson: PROBLEM_SOURCES := $(ERIKSSON)
pure_diffusion_igrm: PROBLEM_SOURCES := $(PURE_DIFFUSION_IGRM)
oil: PROBLEM_SOURCES := $(OIL)
igrm_l2: PROBLEM_SOURCES := $(IGRM_L2)
igrm_heat: PROBLEM_SOURCES := $(IGRM_HEAT)
igrm_eirksson: PROBLEM_SOURCES := $(IGRM_EIRKSSON)

# The problem-local modules have identical names. Always clean and rebuild one
# complete executable before moving to the next problem.
$(PROBLEMS):
	+$(MAKE) --no-print-directory -j1 clean \
		CONFIG="$(CONFIG_PATH)" \
		SOURCE_ALL='$(SOURCES) $(PROBLEM_SOURCES)' EXEC=$@
	+$(MAKE) --no-print-directory -j1 build \
		CONFIG="$(CONFIG_PATH)" \
		SOURCE_ALL='$(SOURCES) $(PROBLEM_SOURCES)' EXEC=$@

problems: $(PROBLEMS)

list-problems:
	@printf '%s\n' $(PROBLEMS)

# Selective cleanup retained for backward compatibility with SOURCE_ALL/EXEC.
clean:
	$(RM) -- $(OBJECT_ALL) $(dep_files) EXEC/$(EXEC)
	$(RM) -- *.o *.mod */*.o */*.mod fort.1

clean-problems:
	$(RM) -- $(addprefix EXEC/,$(PROBLEMS))
	$(RM) -- $(obj_path)/input_data.o $(obj_path)/RHS_fun.o \
		$(obj_path)/igrm_mumps_solver.o $(obj_path)/main.o
	$(RM) -- input_data.mod rhs_fun.mod igrm_mumps_solver.mod

clean-library:
	$(RM) -- $(library_target)

# A library cannot safely reuse Fortran objects produced with another
# compiler, profile, or example configuration. Rebuild its complete object
# set on every public `library` invocation.
clean-library-build:
	$(RM) -r -- $(obj_path)
	$(RM) -- $(library_target)
	$(RM) -- *.o *.mod */*.o */*.mod fort.1

clean-all:
	$(RM) -r -- $(obj_path)
	$(RM) -- $(addprefix EXEC/,$(PROBLEMS)) $(library_target)
	$(RM) -- *.o *.mod */*.o */*.mod fort.1

help:
	@printf '%s\n' \
		'make                         build current SOURCE_ALL/EXEC (default L2)' \
		'make library                 build LIB/libads.a' \
		'make problems                clean-build all eight problem drivers' \
		'make <problem>               clean-build one named driver' \
		'make list-problems' \
		'make clean                   clean current SOURCE_ALL/EXEC selection' \
		'make clean-problems|clean-library|clean-all' \
		'Configuration: $(CONFIG_PATH)'
