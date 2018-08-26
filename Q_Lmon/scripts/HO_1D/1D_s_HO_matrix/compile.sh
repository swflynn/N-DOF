rm *.o, *.mod *.out
gfortran -c sobol_stdnormal.f90
gfortran -c s_sobol.f90
gfortran -O s_sobol.f90 sobol_stdnormal.o
./a.out &
