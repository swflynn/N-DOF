rm *.out *.mod *.o
gfortran -c ndof_mod.f90
gfortran -c ndof.f90
gfortran *.o -o ndof.exe
./ndof.exe
