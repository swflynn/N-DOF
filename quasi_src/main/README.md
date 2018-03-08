# Program Main:
qMC numerical integration of a single water monomer within a water cluster. 

The user defines the number of Basis Functions to describe the monomer with (less than 10) and the highest excitation available to any indivdual basis, and the total excitation for the system. 
The program then computes the Hessian and eigenvectors for the water cluser input geometry. 
A normal mode transformation is then computed and the Potential Difference Matrix (PES - HA) is constructed.
This matrix is then diagonalized (using the LAPACK library) and the fundamental frequencies for the first water within the cluser are computed.

See the development directory for an implementation computing all of the monomers fundamental frequencies. 

## Files:
The following files are necessary for running the program. 

### Makefile
A simple makefile for compiling the program.
This makefile requires both the TIP4P.f90 PES module, and the MBPOL PES module. 
The ssobol.f Fortran code is required for generating the sobol sequences. 
The code is compiled using OMP for parallel implementation. 

### sobol_stdnormal.f90:
This fortran program taken in a vector of sobol points (0,1) and uses the Beasley-Springer-Moro algorithm to transform it to a Gaussian Distribution for numerical integration. 

## run_code
This directory contains an example input file and input geometry for running a simulation.
