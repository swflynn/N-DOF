#=============================================================================80
# The Tip4p and MB-Pol Potential Energy Surfaces are applied for integration.
# See the work of Joel Bowman and Francesco Paesani for the specific surfaces if 
# Interested
# See the MCQMC wiki page for the sobol sequence generator (Fortran; ssobol.f).
#==============================================================================#
OBJS = Tip SSobol SS_std Lmon main
FC = gfortran

all: SSobol Tip SS_std Lmon main

Tip:
	$(FC) -c TIP4P.f90

SSobol:
	$(FC) -c ssobol.f

SS_std:
	$(FC) -c ssobol_stdnormal.f90

Lmon:
	$(FC) -c -fopenmp main.f90 -llapack

main:
	$(FC) -O -Wall -fcheck=all -g -fbacktrace -fopenmp main.o TIP4P.o ssobol.o ssobol_stdnormal.o  mbpol/libmbpol.a   -lstdc++ -llapack

clean:
	rm -rf *.o *.mod *.out $(OBJS)
