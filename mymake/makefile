include m_options
include ./m_files

obj_path = _OBJ
FFLAGS		+= -I./$(obj_path)


EXEC = l2
SOURCE_ALL = $(SOURCES) $(IGRM_L2)
VPATH = $(dir $(SOURCE_ALL))
# $(sort $(dir $(SOURCE_ALL)))

dot_o = .o
objF90 :=  $(notdir $(filter %$(dot_o),$(SOURCE_ALL:.F90=$(dot_o))))
OBJECT_ALL := $(addprefix $(obj_path)/,$(objF90))


all:  EXEC/$(EXEC) 

EXEC/$(EXEC):  $(OBJECT_ALL)
	@$(FF) $(FFLAGS) -o EXEC/$(EXEC) $(OBJECT_ALL) $(USER_LIB)
	@echo DONE



$(obj_path)/%$(dot_o): %.F90
	$(FF) $(FFLAGS) -o $@ -c $<




clean:
	@rm -f $(OBJECT_ALL) $(dep_files) EXEC/$(EXEC)
	rm -f *.o *.mod */*.o */*.mod fort.1


