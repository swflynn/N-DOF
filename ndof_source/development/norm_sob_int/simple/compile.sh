rm *.mod *.o *.out
gfortran -c sobol.f90
gfortran -c sobol_stdnormal.f90
gfortran -c simple.f90
gfortran simple.f90 sobol.o sobol_stdnormal.o
./a.out
