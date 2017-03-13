rm *.mod *.o *.exe
gfortran -c sobol.f90
gfortran -c ndof.f90
gfortran *.o -o ndof.exe
./ndof.exe
