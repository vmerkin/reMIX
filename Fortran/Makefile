F90 = h5fc
#ifort
FFLAGS += -O3 -check bounds -Ddouble_precision -Ddss
#-Dpardiso
#-fast 
#-qopenmp

INCLUDES = -I$(MKLROOT)/include 
LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

FFLAGS += $(INCLUDES) $(LIBS)

.SUFFIXES: .F90 .f .f90

.f90.o:
	$(F90) $(FFLAGS) -c $<

.F90.o:
	$(F90) $(FFLAGS) -c $<

EXE = mix

objects = mixdefs.o\
	mixio.o\
	mixtypes.o\
	mixgeom.o\
	mixsolver.o\
	mix.o

$(EXE): $(objects)
	$(F90) $(FFLAGS) $(LDFLAGS)  -o $(EXE) $(objects) 

clean:
	rm -f *~ *.o *.mod *.optrpt $(EXE)
