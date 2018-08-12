# Program Main:
Code curently takes in the number of spatial dimensions and maximum excitation to evaluate the number of permutations available to the system. 
For each sobol point (standard Sobol Points only), it then evaluates all of the hermite polynomial products and evaluates the potential energy matrix. 
Eigenvalues and matrix elements are provided as a function of iteration to track convergence. 
There is no external field, therefore all integrals are of orthonormal functions.
All diagonal elements should be 1, off-diagonal 0, and the eigenvalues should all be 1. 
Use this test to determine the time and number of sobol points required for highly oscilatory matrix elements with higher excitation and spacial dimensions. 

## Program Files
This current code requires the following files for compilation and execution.

### sobol.f90
Old Fortran code for generating scrambled and non-scrambled sobol points (FSU John Burkardt). 
The code is used to run the `i8_sobol` function to generate a list of sobol points.

### Compile.sh
A sample bash compile file for running the program. Simply run
`$ bash compile.sh` at the terminal to execute.

### uniform_sobol.f90
The main program. 
#### Important Variables:
`m`: Integer, defines the spatial dimension for the integral.

`n`: Integer, number of sobol points to be generated for the entire simulation.
This value must be less than 2**30 for the ssobol.f code to execute. 

`skip`: REAL*8, starting location for the quasi-random sequence.
Suggested to be equal to the value of n.

`start_val`: real, the lower bound for your integral, used for setting the function domain. 

`end_val`: real, the upper bound for your integral, used for setting the function domain. 

`f`: real*8, the integrad evaluations accumulated for MC. 

`r(m,n)`: real, array containing the set of sobol numbers for integration.

`q(m,n)`: real, array containing the set of sobol numbers scaled for your bounds. 
