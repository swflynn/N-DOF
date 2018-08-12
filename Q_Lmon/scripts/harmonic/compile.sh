gfortran -c ssobol.f
gfortran -c *.f90
gfortran -c cg.f
gfortran -O harmonic_approx.o TIP4P.o ssobol.o ssobol_stdnormal.o cg.o 
