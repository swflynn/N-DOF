rm *.mod *.o *.exe
gfortran -c sobol.f90
gfortran -c sobol_stdnormal.f90
gfortran -c int_test.f90
gfortran int_test.f90 sobol.o sobol_stdnormal.o
./a.out
