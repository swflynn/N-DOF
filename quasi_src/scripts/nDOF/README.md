# Program Main:
Code curently takes in the number of spatial dimensions and maximum excitation to evaluate the number of permutations available to the system. 
For each sobol point (standard Sobol Points only), it then evaluates all of the hermite polynomial products and evaluates the potential energy matrix. 
Eigenvalues and matrix elements are provided as a function of iteration to track convergence. 

## Program Files
This current code requires the following files for compilation and execution.

### cg.f
Old Fortran Code used for computing and sorting the eigenvalues associated with the potential energy matrix. 
This code is only used for the `RS` call statement.

### sobol.f90
Old Fortran code for generating scrambled and non-scrambled sobol points. 
The code is used to run the `i8_sobol` function called in the `sobol_stdnormal` subroutine.

### sobol_stdnormal.f90:
This fortran program taken in a vector of sobol points (0,1) and uses the Beasley-Springer-Moro algorithm to transform them to our domain. 
This is used for the `sobol_stdnormal` call statement. 

### Compile.sh
A sample bash compile file for running the program. Simply run
`$ bash compile.sh` at the terminal to execute.

### Ndof.f90
The main program. 
#### Important Variables:
`d`: Integer, defines the spatial dimension, which dictates the length of the sobol point vector. 

`Vmax`: Integer, sets the maximum excitation available to the system for evaluating permutations.
Vmax must be a value between 1,9 or else `permutation` subroutine will not run. 

`Jmax`: Integer containing the total number of permutations, given d and Vmax. 
You program will contain Jmax eigenvalues and a square potential energy matrix of size Jmax.

`v(d,Vmax)`: Integer, contains all of the permutation indices.

`Nsobol`: Integer, number of sobol points to be generated for the entire simulation.
This value must be less than 2**30 for the ssobol.f code to execute. 

`skip`: REAL*8, starting location for the quasi-random sequence.
Suggested to be equal to the value of Nsobol. 

`scrambled_z(d)`: DP, vector containing a unique sobol point for each spatial dimension. 

`herm(Vmax,d)`: DP evaluating each hermite polynomial recursively up to Vmax (index starts from 0). 

`A(Vmax,Vmax,d)`: DP, contains all of the hermite polynomial permutation products. 
This is done for evaluating our integrand as a product over the spatial dimensions. 

`U(Jmax,Jmax)`: DP, Potential Energy Matrix Elements. 

## User Inputs:
The user should choose values for the following:
d, Vmax, Nsobol, and skip.
