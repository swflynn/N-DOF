rm -f *.mod
rm -f *.o
gfortran -c ssobol.f
gfortran -c cg.f
gfortran -c *.f90
gfortran -O quasi_nm.o TIP4P.o ssobol.o ssobol_stdnormal.o cg.o
