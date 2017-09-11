rm *.o *.mod *.out *.dat
gfortran -c sobol_stdnormal.f90
gfortran -c ssobol.f
gfortran -c quas_mc.f90
gfortran -O quasi_mc.f90 ssobol.o sobol_stdnormal.o cg.f -o a.out
