rm *.mod *.o *.exe
gfortran -c sobol.f90
gfortran -c uniform_sobol.f90
gfortran *.o -o unif_sob.exe
./unif_sob.exe
