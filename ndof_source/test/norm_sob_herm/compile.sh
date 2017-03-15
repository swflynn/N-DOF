rm *.mod *.o *.out
gfortran -c sobol.f90
gfortran -c sobol_stdnormal.f90
gfortran -c test.f90
gfortran test.f90 sobol.o sobol_stdnormal.o
./a.out
