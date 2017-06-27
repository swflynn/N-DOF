rm *.o, *.mod *.out
gfortran -c sobol_stdnormal.f90
gfortran -c main.f90
gfortran main.f90 sobol_stdnormal.o
./a.out
