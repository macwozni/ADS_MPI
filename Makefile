include m_files

FC = mpifort
# FFLAGS = -g -pthread -p -fbounds-check
FFLAGS = -pthread -O3

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
	rm -f *.o *.mod */*.o */*.mod

