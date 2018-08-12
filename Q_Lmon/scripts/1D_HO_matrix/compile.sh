rm *.mod *.o *.out *.dat
gfortran -c sobol.f90
gfortran -c sobol_stdnormal.f90
gfortran -c herm_mat.f90
gfortran -O herm_mat.f90 sobol.o sobol_stdnormal.o
