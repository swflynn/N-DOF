rm *.mod *.o *.exe
gfortran -c ndof_mod.f90
gfortran -c ndof.f90
gfortran *.o -o ndof.exe
./ndof.exe
