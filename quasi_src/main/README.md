# Program Main:
qMC numerical integration of a single water monomer within a water cluster. 


The code curently takes in the number of spatial dimensions and maximum excitation to evaluate the number of permutations available to the system. 
For each sobol point (option for non-scrambled, or 3 different scrambling algorithms):
It then evaluates all of the hermite polynomial products and evaluates the potential energy matrix. 
Eigenvalues and matrix elements are provided as a function of iteration to track convergence. 

## Program Files:
This current code requires the following files for compilation and execution.

### quasi_mc.f90:
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

`SAM`: Integer, sets the number of times to perform the scrambled point generation.
This value will repeat your entire calculation, therefore a value of 1 should be used. 
See ssobol.f for details. 

`MAX`: Integer, defining the maximum nmber of digits used for Owen Scrambling. 
Template file suggested setting MAX = 30, see ssobol.f for documentation.

`IFLAG`: Integer, defines the method of scrambling used for generating the Sobol Points. 
Options are 0,1,2,3. 0 = no scrambling, 1 = Owen , 2 = Faure-Tezuka , 3 = Owen + Faure-Tezuka.
From my experience 1, and 3 are the best options for reducing higher-dimensional correlations for numerical integration. 

`scrambled_z(d)`: DP, vector containing a unique sobol point for each spatial dimension. 

`herm(Vmax,d)`: DP evaluating each hermite polynomial recursively up to Vmax (index starts from 0). 

`A(Vmax,Vmax,d)`: DP, contains all of the hermite polynomial permutation products. 
This is done for evaluating our integrand as a product over the spatial dimensions. 

`U(Jmax,Jmax)`: DP, Potential Energy Matrix Elements. 

### cg.f:
Old Fortran Code used for computing and sorting the eigenvalues associated with out potential energy matrix. 
This code is only used for the `RS` call statement.

### ssobol.f:
Old Fortran code for generating scrambled and non-scrambled sobol points. 
The code is used to run the `INSSOBL` and `GOSSOBL` call statements. 

### sobol_stdnormal.f90:
This fortran program taken in a vector of sobol points (0,1) and uses the Beasley-Springer-Moro algorithm to transform them to our domain. 
This is used for the `sobol_stdnormal` call statement. 

### compile.sh:
Bash compile script for running the code as is.
At terminal simply use `$ bash compile.sh` to generate an executable for the script

## User Inputs:
The user should choose values for the following:
d, Vmax, Nsobol, and IFLAG
