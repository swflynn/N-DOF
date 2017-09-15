# Program Main
Code curently takes in the number of spatial dimensions and maximum excitation to evaluate the number of permutations available to the system. 
Code requires a data file containing all of the scrambled sobol points to be used in the calculation. 
Each sobol point is converted to a normal distribution, then evaluated for all of the hermite polynomial products and the potential energy matrix. 
Eigenvalues and matrix elements are provided as a function of iteration to track convergence. 
This code calculated the PE matrix and its eigenvalues for a system without an external field.
Therefore all integrals are of orthonormal functions, off diagonal = 0, on diagonal = 1, eigenvalues = 1. 
Use this code to test methodology before implementing on a real system/ local monomer model. 

## Program Files
This current code requires the following files for compilation and execution.

### s_sobol.m
Matlab script for generating scrambled sobol points and generating the data file `s_sobol_unif.dat`
Define the number of points(Nsobol), the skip (suggested=Nsobol), and the spatial dimension (d).

### cg.f
Old Fortran Code used for computing and sorting the eigenvalues associated with the potential energy matrix. 
This code is only used for the `RS` call statement.

### s_sobol_unif.dat
A data file containing all of the scrambled sobol points to be evaluated (each line contains d points, where d is the spatial dimension. Each point should be seperated by white space. 
This number of sobol points must be specified in the input file.

### sobol_stdnormal.f90
This fortran program taken in a vector of sobol points (0,1) and uses the Beasley-Springer-Moro algorithm to transform them to our domain. 
This is different from the standard sobol method, because we are not generating the sobol points, just changing the distribution. 

### Compile.sh
A sample bash compile file for running the program. Simply run
`$ bash compile.sh` at the terminal to execute.

### Ndof_s.f90
The main program. 

#### Important Variables
`d`: Integer, defines the spatial dimension, which dictates the length of the sobol point vector. 

`Vmax`: Integer, sets the maximum excitation available to the system for evaluating permutations.
Vmax must be a value between 1,9 or else `permutation` subroutine will not run. 

`Jmax`: Integer containing the total number of permutations, given d and Vmax. 
You program will contain Jmax eigenvalues and a square potential energy matrix of size Jmax.

`v(d,Vmax)`: Integer, contains all of the permutation indices.

`Nsobol`: Integer, number of sobol points to be generated for the entire simulation.
This value must be less than 2**30 for the ssobol.f code to execute. 

`scrambled_u(d)`: DP, vector containing a sobol point for each spatial dimension (uniform distribution) 

`scrambled_z(d)`: DP, vector containing a unique sobol point for each spatial dimension (normal distribution).

`herm(Vmax,d)`: DP evaluating each hermite polynomial recursively up to Vmax (index starts from 0). 

`A(Vmax,Vmax,d)`: DP, contains all of the hermite polynomial permutation products. 
This is done for evaluating our integrand as a product over the spatial dimensions. 

`B`: DP, Evaluates product of each hermite polynomial for all the spatial dimensions. 

`U(Jmax,Jmax)`: DP, Potential Energy Matrix Elements. 

## User Inputs:
The user should choose values for the following:
d, Vmax.
Nsobol should be provided in teh input file. 
Also check the convergence analysis to write out at intervals appropriate to your Nsobol. 
