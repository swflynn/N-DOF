rm *.o *.mod *.out *.dat
gfortran -c sobol.f90 
gfortran -c sobol_stdnormal.f90
gfortran -c Ndof.f90
gfortran -O Ndof.f90 sobol.o sobol_stdnormal.o cg.f -o a.out
