rm *.out *.mod *.o
gfortran -c sobol.f90
gfortran -c test.f90
gfortran *.o -o test.exe
./test.exe
