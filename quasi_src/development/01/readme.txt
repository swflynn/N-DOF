gfortran -c *.f90
gfortran -c cg.f
gfortran -O quasi_nm.o TIP4P.o sobol.o sobol_stdnormal.o cg.o mbpol/libmbpol.a   -lstdc++

If you have issue with the make file working this is the compilation order you need
The make file does not work currently just use this method


!=============================================================================================!
 first attempt at the code, this directory can be removed once working code exists
 This code uses the sobol.f90 module
 It is being replaced with ssobol.f, go one directory higher to see this code
 The code in this directory is saved just as a template, it is no longer being developed 9/13/17
!=============================================================================================!
