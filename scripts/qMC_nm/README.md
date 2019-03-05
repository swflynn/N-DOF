# qMC_nm
Code for analyzing a single water monomer (TIP4P PES).
This code uses a normal mode analysis to make the  Hamiltonian Seperable.
It then calculates the Potential Energy Difference matrix and computes the Fundamental Frequencies for the water monomer. 
This calculations is the whole sha-bang! Our First full system calculation. 

The code can output the potential energy difference matrix, the associated eigenvalues and eigenvectors, and the fundamentla frequencies.

# Associated Files

## Compile.sh
A simple bash script for running the code. Simply run the following command on terminal to execute; `$ bash compile.sh`.

## quasi_nm.f90
The main program.

## Input.dat
Line 1: Total Number of Excitations allowed in the system
Line 2: Potential Energy Surface (Only works for TIP4P)
Line 3: Input Water Monomer Geometry, Read in a Hessian Logical
Line 4: Frequency Cutoff (base on the Harmonic Approximation)
Line 5: Freqquency Replace
Line 6: Number of Sobol Points to analze
Line 7: Skip for the sobol point generator
Line 8: Frequency for computing analysis and writing out data to file

Also note: Vmax is hardcoded and needs to be set by the user in the code.
Vmax sets the  highest possible excitation available to any one dimension.
It should contain 9 DOF with excitation values ranging from 0-9 in each DOF.

## ssobol_stdnormal.f90
Fortran program to transform sobol points (0,1) uniformly distributed, to a normal distribution.
We utilize the Beasley-Springer-Moro algorithm to complete this transformation. 

## Running The Code
To run this code you will need an input geometry (xyz) for a water monomer (should be minimized energy).
The code has been simplified to run using the TIP4P potential surface (only), and needs the TIP4P Fortran Module to run.
It also required the Fortran scrambled sobol point generator provided by the MCQMC wiki page. 
Finally we use the cg.f Fortran module for all matrix operations. 
