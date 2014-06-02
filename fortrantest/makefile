FC=gfortran
CFLAGS= -O2

LDIR =

IDIR =
LIBS=-llapack -lblas

PROGRAMS = test01 test02

all: $(PROGRAMS)

test01: nrtype.o nrutil.o nrlinearalg.o test01.o
	$(FC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDIR)

test02: test02.o
	$(FC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDIR)

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
#%: %.o
#	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files

%.o: %.f90
	$(FC) -c $< $(CFLAGS) $(IDIR)

#%: %.o
#	$(FC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDIR)

.PHONY: clean

clean:
	rm -f *.o *.mod *~ core *~ 
	
veryclean: clean
	rm -f *.o *.mod *~ core *~ $(PROGRAMS)
