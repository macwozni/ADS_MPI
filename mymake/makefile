#======================================================================
#
#  The makefile was last modified on $Date: 2006/10/16 18:27:14 $ by $Author: dzubiaur $
#
#======================================================================
#
# To compile:     make C=XXX where XXX is LAPL, LAPLZ, LAPLE (2d code)
#                                         or  R, Z (1d code)
#                            XXX options are set in m_files
#                 make C=XXX OPT=1   - compile in OPT mode
#    preprocess:  make export C=XXX
#    tar:         make tar C=XXX
#    tar & gzip:  make tgz C=XXX
#    clean:       make clean C=XXX <OPT=1> _or_ make realclean
#
#======================================================================


include m_options
  
.SUFFIXES:
.PHONY: build_error banner tar tgz clean realclean link done

SRC_DIR		=  .
REL_DIR		=  ../

# Preprocessing Flags
###############################
include ./m_files
###############################


obj_path = _OBJ
EXEC = l2

# this is real target if nothing is defined
link: banner EXEC/$(EXEC) done


ifeq ($(COMPILER), XLF90)
FFLAGS		+=  -qmoddir=./$(obj_path) -I./$(obj_path)
endif

ifeq ($(COMPILER), IFC)
FFLAGS		+=  -module $(obj_path) -module MODULE -I./$(obj_path) 
endif

ifeq ($(COMPILER), IFORT)
FFLAGS		+=  -module $(obj_path) -module MODULE -I./$(obj_path) 
endif

ifeq ($(COMPILER), PGF90)
FFLAGS		+=  -module $(obj_path) -module MODULE -I./$(obj_path) 
endif

ifeq ($(COMPILER), MPIF90)
FFLAGS		+=  -module $(obj_path) -module MODULE -I./$(obj_path)
endif

ifeq ($(COMPILER), MPICH90)
FFLAGS		+=  -module $(obj_path) -module MODULE -I./$(obj_path)
endif

ifeq ($(COMPILER), GFORTRAN)
FFLAGS		+= -I./$(obj_path) -I./MODULE 
endif

PREC_FLAGS	= $(INCLUDES)
                -DPARALLEL_MODE=$(PARALLEL)

ifeq ($(COMPILER), XLF90)
PREC_FLAGS	= -WF,-DEM_MODE=$(EM),-DC_MODE=$(COMPLEX),-DMAXEQNS_DEF=$(MAXEQNS),-DDIM2=$(DIM2) $(INCLUDES)
endif

# workfiles created by compiler all over the place
TMP_FILES       += *.stb
#======================================================================

#we'll search for source files in these directories
VPATH = $(sort $(dir $(SOURCE_ALL)))

# dot_o can be changed to .obj for Windows
dot_o = .o
objF :=  $(notdir $(filter %$(dot_o),$(SOURCE_ALL:.F90=$(dot_o))))
objc :=  $(notdir $(filter %$(dot_o),$(SOURCE_ALL:.c=$(dot_o))))

OBJECT_ALL := $(addprefix $(obj_path)/,$(objF) $(objc))
#======================================================================
# Link
EXEC/$(EXEC): $(obj_path)/.dummy EXEC/.dummy $(OBJECT_ALL)
	@echo
	@echo Linking EXEC/$(EXEC)
	@echo ----------------------
	@echo $(FF) $(FFLAGS) -o EXEC/$(EXEC) //OBJECT_ALL// $(USER_LIB)
	@$(FF) $(FFLAGS) -o EXEC/$(EXEC) $(OBJECT_ALL) $(USER_LIB)

# initialization banner
banner:
	@echo
ifeq ($(OPT),1)
	@echo OPT compiling into $(obj_path) for $(EXEC)
else
	@echo Compiling into $(obj_path) for $(EXEC)
endif
	@echo

done:
	@echo
	@echo Done: EXEC/$(EXEC) complete.
	@echo

$(obj_path)/.dummy:
	@echo Creating $(obj_path) directory ...
	@if [ -d  $(obj_path) ] ; then \
		touch $@; \
	   else mkdir $(obj_path); touch $@ ; \
	   fi

EXEC/.dummy:
	@echo Creating EXEC directory ...
	@if [ -d  EXEC ] ; then \
		touch $@; \
	   else mkdir EXEC; touch $@ ; \
	   fi


# Compile .f90 to create .F90
$(obj_path)/%$(dot_o): %.F90
	@echo Compiling $<
	@echo $(FF) $(FFLAGS) $(PREC_FLAGS) $(INCLUDES) -o $@ -c $<
	$(FF) $(FFLAGS) $(PREC_FLAGS) $(INCLUDES) -o $@ -c $<

$(obj_path)/%$(dot_o): %.c
	@echo Compiling $<
	$(CC) -c $< $(CFLAGS) -o $@


# create default compiler options from the template
m_options: MAKE-SHARED/m_compile
	@if [ -f m_options ] ; then \
	    touch m_options ; \
	    echo "MAKE-SHARED/m_compile changed !" ;\
	    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ !" ;\
	  else \
	    cp MAKE-SHARED/m_compile m_options ;\
	  fi
#==================================================================
# auto dependencies magic
ifndef NODEPENDECIES

# auto dependencies require gcc, if you do not have it, then you
#  can desable autodependenceis using
#     make NODEPENDECIES=1 ...other options....

dep_files = $(patsubst %$(dot_o),%.P,$(OBJECT_ALL))

-include $(dep_files) $(obj_path)/.dummy
$(dep_files): $(obj_path)/.dummy

FIX_P = ( echo -n $@ $(obj_path)/ ; cat $@-tmp ) > $@

$(obj_path)/%.P: %.c
	@echo mkdeps $< '>' $@ 
	@(gcc -MM $(PREC_FLAGS) $(INCLUDES) $< > $@-tmp) 2> /dev/null 
	@$(FIX_P)
	@rm $@-tmp

$(obj_path)/%.P: %.F90
	@echo mkdeps $< '>' $@ 
	@sed \
	    -e '/^[^#]/d' \
	    -e '/.*\\ *$$/s/\\* *$$//' \
	    -e '/\/\*/s/\/\*.*//' \
	    -e "s/'//g" \
	    -e 's/`//g' \
		  		$<  > $(*F)-sed.c
	@(gcc -xc -MM $(PREC_FLAGS) $(INCLUDES) $(*F)-sed.c > $@-tmp) 
	@echo $@ " " \\ > $@
	@sed \
	      -e 's,$(*F)-sed\.o,$$\(obj_path\)/$(*F)$(dot_o),g' \
	      -e 's,$(*F)-sed.c,$(*F).F90,g' \
	      -e 's,$(*F)-sed:,$$\(obj_path\)/$(*F)$(dot_o):,g' \
	      -e 's,$(*F)-sed$$,$(*F).F90,g' \
	   $@-tmp >> $@
	@rm -f $(*F)-sed.c $@-tmp

endif
#======================================================================
#export magic
OTHER_FILES += makefile m_files \
   MAKE-SHARED/makefile-common MAKE-SHARED/m_compile MAKE-SHARED/prep.tcl
exp_path = ./export_$(C)
blk_path = $(sort $(dir $(BLK_FILES)))
# for ffld/commons we need to add also ffld
blk_exp = $(blk_path) $(dir $(patsubst %/,%,$(blk_path)))
all_dirs = $(patsubst ./%,%,$(dir $(SOURCE_ALL) $(OTHER_FILES)) $(blk_exp))
exp_dirs = $(sort $(addprefix $(exp_path)/,$(all_dirs)))
export:
	rm -rf $(exp_path)
	mkdir $(exp_path)
	for d in $(exp_dirs) ; do mkdir $$d ; done
	for f in $(SOURCE_ALL) $(BLK_FILES) ; do echo $(exp_path)/$$f ;\
	  (tclsh MAKE-SHARED/prep.tcl $(COMPLEX) $(EM) $(MAXEQNS) \
                            $(PARALLEL)\
                            < $$f > $(exp_path)/$$f); \
	  done
	for f in $(OTHER_FILES) ; do echo $(exp_path)/$$f ;\
	  cp $$f $(exp_path)/$$f ; done
	@echo ...
	tar -czf $(TAR_NAME)_$(C).tgz $(exp_path)
	@echo $(exp_path)/ AND $(TAR_NAME)_$(C).tgz CREATEDLIBIRC

# Ftof_FILES	:= $(patsubst %.F90,%.f90,$(wildcard */*.F90))

# tar or tar/gzip

tar: 
	@tar -cvf $(TAR_NAME).tar $(TAR_FILES) $(OTHER_FILES)
tgz: 
	@tar -czvf $(TAR_NAME).tgz $(TAR_FILES) $(OTHER_FILES) 

list:
	@echo LIST OF SOURCE FILES:
	@echo $(SOURCE_ALL)
	@echo
	@echo LIST OF SOURCE PATHS SEARCHED IN:
	@echo $(VPATH)
	@echo
	@echo LIST OF DIRECTORIES SEARCHED FOR INCLUDES:
	@echo $(INCLUDES)
	@echo
	@echo PREPROCESSOR AND COMPILER OPTIONS:
	@echo $(PREC_FLAGS)
	@echo
	@echo $(FFLAGS)
	@echo

clean:
	@rm -f $(OBJECT_ALL) $(dep_files) EXEC/$(EXEC)
	rm -f *.o *.mod */*.o */*.mod fort.1

#realclean:
#	@rm -rf _OBJ_* EXEC export* $(TAR_NAME)* $(TMP_FILES)
#   these are hidden files created by editors and cvs conflicts
#	@rm -f .\#* */.\#* */*/.\#* \#*\# */\#*\# */*/\#*\# *~ */*~ */*/*~

