rm *.mod *.o *.out
gfortran -c sobol.f90
gfortran -c unif_herm.f90
gfortran *.o -o unif.out
./unif.out
