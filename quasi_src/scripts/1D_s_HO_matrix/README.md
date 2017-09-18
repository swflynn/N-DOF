# 1D_s_HO_matrix
Code for evaluating the Potential Energy Matrix (QMC Scrambled-Sobol Sequence) for a Harmonic Oscillator Wavefunctions in a single spatial dimension.
Code takes in what degree Hermite Polynomial you would like to calculate up to (first 10 available), and the number of Sobol Points to use for the evaluation (each sobol point is converted from a uniform distribution to a normal distribution). 
The points must be generated in a seperate program, and must be given as input. 
Then calculate the PE matrix elements for each permutation of two wavefunctions up to deg. 
This code represents a single basis function to be used in the n-Dimensional code (i.e. a basis is composed of a product of 2 wavefunctions and a gaussian), a single spatial dimension. 
This code calculated the PE matrix for a system without an external field.
Therefore all integrals are orthonormal, off diagonal = 0, on diagonal = 1.

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


`norm(d, Nsobol)`: DP, vector containing normal distribution sobol point. 

`coef`: DP, Coefficient to multiply Hermite Polynomial by to get Wavefunction. 

`herm`: DP, evaluating each hermite polynomial recursively up to deg.

`A(deg,deg)`: DP, Potential Energy Matrix Elements for a single dimension (matrix is orthonormal, diagonal=1, off-diagonal=0). 
