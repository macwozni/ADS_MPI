include m_options
include ./m_files

.DEFAULT_GOAL := all
.PHONY: all clean

obj_path = _OBJ
FFLAGS		+= -J. -I. -I./$(obj_path)


EXEC = l2
SOURCE_ALL = $(SOURCES) $(L2)
VPATH = $(dir $(SOURCE_ALL))
# $(sort $(dir $(SOURCE_ALL)))

dot_o = .o
objF90 :=  $(notdir $(filter %$(dot_o),$(SOURCE_ALL:.F90=$(dot_o))))
OBJECT_ALL := $(addprefix $(obj_path)/,$(objF90))


all:  EXEC/$(EXEC)

EXEC/$(EXEC):  $(OBJECT_ALL) | EXEC
	@$(FF) $(FFLAGS) -o EXEC/$(EXEC) $(OBJECT_ALL) $(USER_LIB)
	@echo DONE



$(obj_path)/%$(dot_o): %.F90 | $(obj_path)
	$(FF) $(FFLAGS) -o $@ -c $<



EXEC $(obj_path):
	mkdir -p $@




clean:
	@rm -f $(OBJECT_ALL) $(dep_files) EXEC/$(EXEC)
	rm -f *.o *.mod */*.o */*.mod fort.1
