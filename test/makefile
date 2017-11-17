FORTRAN=ifort
FLAGS=-O0

SRCDIR=../src
TESTDIR=.
FRUITDIR=../fruit/src

OBJS=knot_vector.o
TESTS=knot_vector_test.o fruit_driver.o
FRUITS=fruit_util.o fruit.o

.PHONY : all
all : test


tests : $(FRUITS) $(OBJS) $(TESTS)
	$(FORTRAN) $(FLAGS) $^ -o $@

%.o : $(SRCDIR)/%.F90
	$(FORTRAN) $(FLAGS) -c $< -o $@
%.o : $(TESTDIR)/%.F90
	$(FORTRAN) $(FLAGS) -c $< -o $@
%.o : $(TESTDIR)/%.F90
	$(FORTRAN) $(FLAGS) -c $< -o $@
%.o : $(FRUITDIR)/%.F90
	$(FORTRAN) $(FLAGS) -c $< -o $@
	

.PHONY : test
test : tests
	./$<

.PHONY : clean
clean :
	rm -f tests
	rm -f result*.xml
	rm -f *.mod 
	rm -f *.o
	rm -f *~
	rm -f src/*~
	rm -f test/*~
