rm *.mod *.o *.exe
gfortran -c Ndof_s_module.f90
gfortran -c Ndof_s.f90
gfortran *.o -o Ndof.exe
./Ndof.exe
