CXX = g++									# program for compiling c++, g++ = default
FC = gfortran							# program for compiling Fortran programs, default = F77
CXXFLAGS = -O2 -Wall			# extra flags for c++ compiler
FFLAGS = -O2 -Wall				# extra flags for fortran compiler

LD = $(CXX)								# Libraries needed for linking (the linker)
LDFLAGS = $(CXXFLAGS)			# extra flags for linker

all: quasi_nm							# first target run with make,

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^

HD_nDOF: quasi_nm.o TIP4P.o cg.o sobol.o sobol_stdnormal.o \
	mbpol/libmbpol.a
	$(FC) -o $@ $^ mbpol/libmbpol.a $(FFLAGS) -lstdc++

clean:
	rm -f libX2O.a *.o *.mod quasi_nm
