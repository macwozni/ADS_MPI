# Isolated adapter for the historical SOURCE_ALL/EXEC interface.
#
# This file deliberately owns no official source list.  It compiles only the
# ordered paths explicitly supplied by the caller and keeps every intermediate
# file below BUILD_ROOT/_LEGACY_OBJ/EXEC.

ifndef BUILD_ROOT_PATH
$(error legacy-source-build.mk requires BUILD_ROOT_PATH)
endif

EXEC ?= l2

LEGACY_SOURCE_FILES := $(strip $(SOURCE_ALL))
LEGACY_SOURCE_PATHS := $(foreach source,$(LEGACY_SOURCE_FILES),$(abspath $(source)))
LEGACY_SOURCE_BASENAMES := $(notdir $(LEGACY_SOURCE_PATHS))
LEGACY_DUPLICATE_BASENAMES := $(sort $(foreach name,$(LEGACY_SOURCE_BASENAMES),\
	$(if $(word 2,$(filter $(name),$(LEGACY_SOURCE_BASENAMES))),$(name))))
LEGACY_UNSUPPORTED_SOURCES := $(filter-out %.F90,$(LEGACY_SOURCE_PATHS))

LEGACY_ROOT_DIR := $(BUILD_ROOT_PATH)/_LEGACY_OBJ
LEGACY_OBJ_DIR := $(LEGACY_ROOT_DIR)/$(EXEC)
LEGACY_EXEC_DIR := $(BUILD_ROOT_PATH)/EXEC
LEGACY_EXECUTABLE := $(LEGACY_EXEC_DIR)/$(EXEC)
LEGACY_ROOT_MARKER := $(LEGACY_ROOT_DIR)/.ads-legacy-root
LEGACY_OBJ_MARKER := $(LEGACY_OBJ_DIR)/.ads-legacy-object-dir
LEGACY_EXEC_MARKER := $(LEGACY_EXEC_DIR)/.ads-legacy-exec-dir
LEGACY_BUILD_STAMP := $(LEGACY_OBJ_DIR)/.ads-legacy-build
LEGACY_OBJECTS := $(addprefix $(LEGACY_OBJ_DIR)/,\
	$(addsuffix .o,$(basename $(LEGACY_SOURCE_BASENAMES))))

MODULE_OUTPUT ?= -J.
ifneq ($(filter -J.%,$(MODULE_OUTPUT)),)
LEGACY_MODULE_FLAGS := -J$(LEGACY_OBJ_DIR)
else ifneq ($(filter -module,$(MODULE_OUTPUT)),)
LEGACY_MODULE_FLAGS := -module $(LEGACY_OBJ_DIR)
else
LEGACY_MODULE_FLAGS := $(MODULE_OUTPUT)
endif

LEGACY_COMPILE_FLAGS := $(FFLAGS) $(LEGACY_MODULE_FLAGS) \
	-I$(LEGACY_OBJ_DIR) -I$(BUILD_ROOT_PATH)

.PHONY: legacy-build legacy-clean legacy-clean-all legacy-force \
	legacy-prepare-exec-dir \
	validate-legacy-paths validate-legacy-input

legacy-force:

validate-legacy-paths:
	@case '$(EXEC)' in \
		''|[-.]*|*[!A-Za-z0-9_.-]*) \
			printf 'Unsafe legacy EXEC: %s\nUse letters, digits, _, . and -; start with a letter, digit, or _.\n' \
				'$(EXEC)' >&2; exit 2 ;; \
	esac
	@case '$(BUILD_ROOT_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(SRC_DIR)'|'$(PROBLEMS_DIR)') \
			printf 'Unsafe legacy BUILD_ROOT: %s\n' '$(BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac
	@if test -L '$(LEGACY_ROOT_DIR)' || test -L '$(LEGACY_OBJ_DIR)' || \
		test -L '$(LEGACY_EXEC_DIR)'; then \
		printf 'Refusing symlinked legacy object path: %s\n' \
			'$(LEGACY_OBJ_DIR)' >&2; exit 2; \
	fi

validate-legacy-input: validate-legacy-paths
	@if test -z '$(LEGACY_SOURCE_FILES)'; then \
		printf '%s\n' 'SOURCE_ALL was supplied but is empty.' >&2; exit 2; \
	fi
	@if test -n '$(LEGACY_UNSUPPORTED_SOURCES)'; then \
		printf 'SOURCE_ALL supports ordered .F90 files only; unsupported: %s\n' \
			'$(LEGACY_UNSUPPORTED_SOURCES)' >&2; exit 2; \
	fi
	@if test -n '$(LEGACY_DUPLICATE_BASENAMES)'; then \
		printf 'SOURCE_ALL contains duplicate basenames: %s\n' \
			'$(LEGACY_DUPLICATE_BASENAMES)' >&2; exit 2; \
	fi

$(LEGACY_ROOT_MARKER): validate-legacy-paths
	@legacy_root='$(LEGACY_ROOT_DIR)'; marker='$@'; \
	if test -e "$$legacy_root" && test ! -d "$$legacy_root"; then \
		printf 'Legacy object root is not a directory: %s\n' "$$legacy_root" >&2; exit 2; \
	fi; \
	if test -e "$$marker"; then \
		if test ! -f "$$marker" || test -L "$$marker" || \
			! grep -Fqx -- 'ADS_MPI_LEGACY_ROOT=$(LEGACY_ROOT_DIR)' "$$marker"; then \
			printf 'Refusing invalid legacy root marker: %s\n' "$$marker" >&2; exit 2; \
		fi; \
		exit 0; \
	fi; \
	if test -d "$$legacy_root" && test ! -e "$$marker" && \
		test -n "$$(find "$$legacy_root" -mindepth 1 -maxdepth 1 -print -quit)"; then \
		printf 'Refusing to adopt non-empty unmarked directory: %s\n' "$$legacy_root" >&2; exit 2; \
	fi; \
	mkdir -p "$$legacy_root"; \
	marker_tmp="$$marker.tmp"; \
	printf 'ADS_MPI_LEGACY_ROOT=%s\n' "$$legacy_root" > "$$marker_tmp"; \
	mv -f -- "$$marker_tmp" "$$marker"

$(LEGACY_OBJ_MARKER): $(LEGACY_ROOT_MARKER) validate-legacy-paths
	@legacy_dir='$(LEGACY_OBJ_DIR)'; marker='$@'; \
	if test -e "$$legacy_dir" && test ! -d "$$legacy_dir"; then \
		printf 'Legacy object path is not a directory: %s\n' "$$legacy_dir" >&2; exit 2; \
	fi; \
	if test -e "$$marker"; then \
		if test ! -f "$$marker" || test -L "$$marker" || \
			! grep -Fqx -- 'ADS_MPI_LEGACY_EXEC=$(EXEC)' "$$marker" || \
			! grep -Fqx -- 'ADS_MPI_LEGACY_DIR=$(LEGACY_OBJ_DIR)' "$$marker"; then \
			printf 'Refusing invalid legacy object marker: %s\n' "$$marker" >&2; exit 2; \
		fi; \
		exit 0; \
	fi; \
	if test -d "$$legacy_dir" && test ! -e "$$marker" && \
		test -n "$$(find "$$legacy_dir" -mindepth 1 -maxdepth 1 -print -quit)"; then \
		printf 'Refusing to adopt non-empty unmarked directory: %s\n' "$$legacy_dir" >&2; exit 2; \
	fi; \
	mkdir -p "$$legacy_dir"; \
	marker_tmp="$$marker.tmp"; \
	{ \
		printf 'ADS_MPI_LEGACY_EXEC=%s\n' '$(EXEC)'; \
		printf 'ADS_MPI_LEGACY_DIR=%s\n' "$$legacy_dir"; \
	} > "$$marker_tmp"; \
	mv -f -- "$$marker_tmp" "$$marker"

legacy-prepare-exec-dir: validate-legacy-paths
	@exec_dir='$(LEGACY_EXEC_DIR)'; marker='$(LEGACY_EXEC_MARKER)'; \
	if test -e "$$exec_dir" && test ! -d "$$exec_dir"; then \
		printf 'Legacy executable path is not a directory: %s\n' "$$exec_dir" >&2; exit 2; \
	fi; \
	if test -e "$$marker"; then \
		if test ! -f "$$marker" || test -L "$$marker" || \
			! grep -Fqx -- 'ADS_MPI_LEGACY_EXEC_DIR=$(LEGACY_EXEC_DIR)' "$$marker"; then \
			printf 'Refusing invalid legacy executable marker: %s\n' "$$marker" >&2; exit 2; \
		fi; \
		exit 0; \
	fi; \
	if test -d "$$exec_dir"; then exit 0; fi; \
	mkdir -p "$$exec_dir"; \
	marker_tmp="$$marker.tmp"; \
	printf 'ADS_MPI_LEGACY_EXEC_DIR=%s\n' "$$exec_dir" > "$$marker_tmp"; \
	mv -f -- "$$marker_tmp" "$$marker"

$(LEGACY_BUILD_STAMP): legacy-force | $(LEGACY_OBJ_MARKER)
	@stamp_tmp='$@.tmp'; \
	{ \
		printf '%s\n' 'CONFIG=$(CONFIG_PATH)'; \
		printf '%s\n' 'EXEC=$(EXEC)'; \
		printf '%s\n' 'SOURCES=$(LEGACY_SOURCE_PATHS)'; \
		printf '%s\n' 'FF=$(FF)'; \
		printf '%s\n' 'FFLAGS=$(LEGACY_COMPILE_FLAGS)'; \
		printf '%s\n' 'USER_LIB=$(USER_LIB)'; \
	} > "$$stamp_tmp"; \
	if test -r '$@' && cmp -s "$$stamp_tmp" '$@'; then \
		$(RM) -- "$$stamp_tmp"; \
	else \
		mv -f -- "$$stamp_tmp" '$@'; \
	fi

# Emit one explicit compile rule per caller-supplied source.  Rules are only
# generated for a structurally valid list, so duplicate object targets cannot
# silently override one another before validation reports the error.
ifneq ($(LEGACY_SOURCE_FILES),)
ifeq ($(LEGACY_DUPLICATE_BASENAMES),)
ifeq ($(LEGACY_UNSUPPORTED_SOURCES),)
define DEFINE_LEGACY_COMPILE_RULE
$(LEGACY_OBJ_DIR)/$(basename $(notdir $(1))).o: $(1) $(MAKEFILE_LIST) | $(LEGACY_OBJ_MARKER)
	$$(FF) $$(LEGACY_COMPILE_FLAGS) -o "$$@" -c "$$<"
endef
$(foreach source,$(LEGACY_SOURCE_PATHS),$(eval $(call DEFINE_LEGACY_COMPILE_RULE,$(source))))

$(LEGACY_OBJECTS): $(LEGACY_BUILD_STAMP)

# Materialize the user's order as dependencies for Fortran module consumers.
define add-legacy-ordered-dependencies
$(if $(word 2,$1),$(eval $(word 2,$1): $(word 1,$1))$(call add-legacy-ordered-dependencies,$(wordlist 2,$(words $1),$1)))
endef
$(call add-legacy-ordered-dependencies,$(LEGACY_OBJECTS))

$(LEGACY_EXECUTABLE): $(LEGACY_OBJECTS) $(LEGACY_BUILD_STAMP) | legacy-prepare-exec-dir
	$(FF) $(LEGACY_COMPILE_FLAGS) -o "$@" $(LEGACY_OBJECTS) $(USER_LIB)
	@printf 'Built legacy executable %s\n' '$@'

legacy-build: validate-legacy-input $(LEGACY_EXECUTABLE)
else
legacy-build: validate-legacy-input
endif
else
legacy-build: validate-legacy-input
endif
else
legacy-build: validate-legacy-input
endif

legacy-clean: validate-legacy-input
	@legacy_dir='$(LEGACY_OBJ_DIR)'; marker='$(LEGACY_OBJ_MARKER)'; \
	exec_dir='$(LEGACY_EXEC_DIR)'; exec_marker='$(LEGACY_EXEC_MARKER)'; \
	legacy_owned=0; \
	if test -e "$$exec_marker" && \
		{ test ! -f "$$exec_marker" || test -L "$$exec_marker" || \
		! grep -Fqx -- 'ADS_MPI_LEGACY_EXEC_DIR=$(LEGACY_EXEC_DIR)' "$$exec_marker"; }; then \
		printf 'Refusing invalid legacy executable marker: %s\n' \
			"$$exec_marker" >&2; exit 2; \
	fi; \
	if test -e "$$legacy_dir"; then \
		if test ! -d "$$legacy_dir" || test -L "$$legacy_dir" || \
			test ! -f "$$marker" || test -L "$$marker" || \
			! grep -Fqx -- 'ADS_MPI_LEGACY_EXEC=$(EXEC)' "$$marker" || \
			! grep -Fqx -- 'ADS_MPI_LEGACY_DIR=$(LEGACY_OBJ_DIR)' "$$marker"; then \
			printf 'Refusing to remove unmarked legacy directory: %s\n' \
				"$$legacy_dir" >&2; exit 2; \
		fi; \
		legacy_owned=1; \
		$(RM) -r -- "$$legacy_dir"; \
	fi; \
	if test "$$legacy_owned" -eq 1; then \
		$(RM) -- '$(LEGACY_EXECUTABLE)'; \
	fi; \
	if test -f "$$exec_marker" && \
		test -z "$$(find "$$exec_dir" -mindepth 1 -maxdepth 1 \
			! -name '.ads-legacy-exec-dir' -print -quit)"; then \
		$(RM) -- "$$exec_marker"; \
		rmdir -- "$$exec_dir"; \
	fi

# Validate every registered child before deleting anything.  This makes a
# corrupted or foreign directory fail closed instead of being partly removed.
legacy-clean-all:
	@legacy_root='$(LEGACY_ROOT_DIR)'; root_marker='$(LEGACY_ROOT_MARKER)'; \
	exec_dir='$(LEGACY_EXEC_DIR)'; exec_marker='$(LEGACY_EXEC_MARKER)'; \
	if test ! -e "$$legacy_root"; then exit 0; fi; \
	case '$(BUILD_ROOT_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(SRC_DIR)'|'$(PROBLEMS_DIR)') \
			printf 'Unsafe legacy BUILD_ROOT: %s\n' '$(BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac; \
	if test ! -d "$$legacy_root" || \
		test -L "$$legacy_root" || test ! -f "$$root_marker" || \
		test -L "$$root_marker" || \
		! grep -Fqx -- 'ADS_MPI_LEGACY_ROOT=$(LEGACY_ROOT_DIR)' "$$root_marker"; then \
		printf 'Refusing to remove unmarked legacy root: %s\n' "$$legacy_root" >&2; exit 2; \
	fi; \
	if test -L "$$exec_dir" || \
		{ test -e "$$exec_marker" && \
		  { test ! -f "$$exec_marker" || test -L "$$exec_marker" || \
		    ! grep -Fqx -- 'ADS_MPI_LEGACY_EXEC_DIR=$(LEGACY_EXEC_DIR)' "$$exec_marker"; }; }; then \
		printf 'Refusing invalid legacy executable directory: %s\n' \
			"$$exec_dir" >&2; exit 2; \
	fi; \
	if test -n "$$(find "$$legacy_root" -mindepth 1 -maxdepth 1 \
		! -name '.ads-legacy-root' ! -type d -print -quit)" || \
		test -n "$$(find "$$legacy_root" -mindepth 1 -maxdepth 1 \
		-type d -name '.*' -print -quit)"; then \
		printf 'Refusing unexpected entries in legacy root: %s\n' "$$legacy_root" >&2; exit 2; \
	fi; \
	for legacy_dir in "$$legacy_root"/*; do \
		test -e "$$legacy_dir" || continue; \
		exec_name=$${legacy_dir##*/}; marker="$$legacy_dir/.ads-legacy-object-dir"; \
		case "$$exec_name" in ''|[-.]*|*[!A-Za-z0-9_.-]*) \
			printf 'Unsafe registered legacy EXEC: %s\n' "$$exec_name" >&2; exit 2 ;; esac; \
		if test ! -d "$$legacy_dir" || test -L "$$legacy_dir" || \
			test ! -f "$$marker" || test -L "$$marker" || \
			! grep -Fqx -- "ADS_MPI_LEGACY_EXEC=$$exec_name" "$$marker" || \
			! grep -Fqx -- "ADS_MPI_LEGACY_DIR=$$legacy_dir" "$$marker"; then \
			printf 'Refusing unmarked legacy child: %s\n' "$$legacy_dir" >&2; exit 2; \
		fi; \
	done; \
	for legacy_dir in "$$legacy_root"/*; do \
		test -e "$$legacy_dir" || continue; \
		exec_name=$${legacy_dir##*/}; \
		$(RM) -- '$(LEGACY_EXEC_DIR)'/"$$exec_name"; \
	done; \
	$(RM) -r -- "$$legacy_root"; \
	if test -f "$$exec_marker" && \
		test -z "$$(find "$$exec_dir" -mindepth 1 -maxdepth 1 \
			! -name '.ads-legacy-exec-dir' -print -quit)"; then \
		$(RM) -- "$$exec_marker"; \
		rmdir -- "$$exec_dir"; \
	fi
