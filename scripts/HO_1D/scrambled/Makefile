OBJS = ss_std HO main
FC = gfortran

all: ss_std HO main

ss_std:
	$(FC) -c ss_stdnorm.f90

HO:
	$(FC) -c HO_1D_s.f90

main:
	$(FC) -O HO_1D_s.o ss_stdnorm.o

clean:
	rm -rf *.o *.mod *.out $(OBJS)
