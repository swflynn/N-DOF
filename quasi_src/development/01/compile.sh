rm a.out
rm matrix.dat
rm eigenvalues.dat
rm output.dat
gfortran -c *.f90
gfortran -c cg.f
gfortran -O quasi_nm.o TIP4P.o sobol.o sobol_stdnormal.o cg.o mbpol/libmbpol.a   -lstdc++
