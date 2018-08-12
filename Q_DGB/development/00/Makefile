OBJS = Sobol S_std DGB main
FC = gfortran

all: Sobol S_std DGB main

Sobol:
	$(FC) -c sobol.f90

S_std:
	$(FC) -c sobol_stdnormal.f90

DGB:
	$(FC) -c DGB.f90 -llapack

main:
	$(FC) -O DGB.o sobol.o sobol_stdnormal.o -llapack

clean:
	rm -rf *.o *.mod *.out $(OBJS)
