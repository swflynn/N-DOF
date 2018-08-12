# 1D_HO_matrix
Code for evaluating the Potential Energy Matrix (QMC Sobol Sequence) for a Harmonic Oscillator Wavefunctions in a single spatial dimension.
Code takes in what degree Hermite Polynomial you would like to calculate up to (first 10 available), and the number of Sobol Points to use for the evaluation (each sobol point is converted from a uniform distribution to a normal distribution). 
Then calculate the PE matrix elements for each permutation of two wavefunctions up to deg. 
This code represents a single basis function to be used in the n-Dimensional code (i.e. a basis is composed of a product of 2 wavefunctions and a gaussian), a single spatial dimension. 
This code calculated the PE matrix for a system without an external field.
Therefore all integrals are orthonormal, off diagonal = 0, on diagonal = 1.

## Associated Files

### herm_mat.f90
The program main

### sobol_stdnormal.f90
This fortran program taken in a vector of sobol points (0,1) uniformly distributed, and uses the Beasley-Springer-Moro algorithm to transform them to a normal distribution. 

### sobol.f90
Fortran code for generating scrambled and non-scrambled sobol points (see FSU John Burkardt for src code). 
The code is used to run the i8_sobol function called in the sobol_stdnormal subroutine.

### Compile.sh
A sample bash compile file for running the program. Simply run
`$ bash compile.sh` at the terminal to execute.
