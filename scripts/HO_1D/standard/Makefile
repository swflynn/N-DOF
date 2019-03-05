OBJS = Sobol S_std HO main
FC = gfortran

all: Sobol S_std HO main

Sobol:
	$(FC) -c sobol.f90

S_std:
	$(FC) -c sobol_stdnormal.f90

HO:
	$(FC) -c HO_1D.f90

main:
	$(FC) -O HO_1D.o sobol.o sobol_stdnormal.o

clean:
	rm -rf *.o *.mod *.out $(OBJS)
