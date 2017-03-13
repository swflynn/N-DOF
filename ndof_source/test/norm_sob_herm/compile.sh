rm *.mod *.o *.out
gfortran -c my_lib.f90
gfortran -c main.f90
gfortran *.o -o main_prog.out
./main_prog.out
