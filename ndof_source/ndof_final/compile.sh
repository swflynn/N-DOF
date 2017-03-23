rm *.mod *.o *.out *.dat
gfortran -c sobol.f90
gfortran -c sobol_stdnormal.f90
gfortran -c ndof.f90
gfortran ndof.f90 sobol.o sobol_stdnormal.o
./a.out
