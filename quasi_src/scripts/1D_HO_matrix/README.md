# 1D_HO_matrix
Code for evaluating the Potential Energy Matrix (QMC Sobol Sequence) for a Harmonic Oscillator Wavefunctions in a single spatial dimension.
Code takes in what degree Hermite Polynomial you would like to calculate up to (first 10 available), and the number of Sobol Points to use for the evaluation (each sobol point is converted from a uniform distribution to a normal distribution). 
It then calculates the matrix elements for each permutation of two wavefunctions up to deg. 
The code represents a single basis function to be used in the n-Dimensional code (i.e. a basis is composed of a product of 2 wavefunctions and a gaussian). 
This code calculated the PE matrix for a system without an external field.
Therefore all integrals are of orthonormal functions, off diagonal = 0, on diagonal = 1, eigenvalues = 1. 

## Program Files
This current code requires the following files for compilation and execution.

### sobol_stdnormal.f90
This fortran program taken in a vector of sobol points (0,1) uniformly distributed, and uses the Beasley-Springer-Moro algorithm to transform them to a normal distribution. 

### sobol.f90
Old Fortran code for generating scrambled and non-scrambled sobol points (FSU John Burkardt). 
The code is used to run the i8_sobol function called in the sobol_stdnormal subroutine.

### Compile.sh
A sample bash compile file for running the program. Simply run
`$ bash compile.sh` at the terminal to execute.

### herm_mat.f90
The main program.

#### Important Variables
`Nsobol`: Integer, number of sobol points to be generated for the entire simulation.

`deg`: Integer, defines the order of Hermite Polynomial to calculate up to (index starts at 1 not 0), valid for the first ten polynomials. 

`scrambled_u(d)`: DP, vector containing a sobol point for each spatial dimension (uniform distribution) 

`scrambled_z(d)`: DP, vector containing a unique sobol point for each spatial dimension (normal distribution).

`coef`: DP, Coefficient to multiply Hermite Polynomial by to get Wavefunction. 

`herm`: DP, evaluating each hermite polynomial recursively up to deg.

`A(deg,deg)`: DP, contains all of the hermite polynomial permutation products. 
This is done for evaluating our integrand as a product over the spatial dimensions. 

`B`: DP, Evaluates product of each hermite polynomial for all the spatial dimensions. 

`U(Jmax,Jmax)`: DP, Potential Energy Matrix Elements. 

## User Inputs:
The user should choose values for the following:
Nsobol, skip (suggested to be equal to Nsobol), and deg.
