include m_files

FC = /opt/mpich/bin/mpif90
#FC = /opt/mpich/bin/mpif90 -f90=gfortran
# FFLAGS = -g -pthread -p -fbounds-check
FFLAGS = -pthread -g
#-DIPRINT
#-DIDEBUG

OS := $(SOURCES:.F90=.o)
OBJECTS := $(OS:.f=.o)
OBJECTS_BLAS := $(SOURCES_BLAS:.f=.o)
OBJECTS_LAPACK := $(SOURCES_LAPACK:.f=.o)

SOURCES += $(SOURCES_LAPACK)
OBJECTS += $(OBJECTS_LAPACK)
SOURCES += $(SOURCES_BLAS)
OBJECTS += $(OBJECTS_BLAS)

l2: $(OBJECTS) 
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^

.PHONY: clean
clean:
	rm -f *.o *.mod */*.o */*.mod l2

