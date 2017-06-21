rm *.o, *.mod *.out
gfortran -c n_s_stdnorm.f90
gfortran -c n_s_sobol.f90
gfortran n_s_sobol.f90 n_s_stdnorm.o
./a.out
rm *.o *.mod a.out
