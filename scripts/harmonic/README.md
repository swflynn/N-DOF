# Harmonic
This code computes the Fundamental Frequencies for a water cluster using the TIP4P potential only.
It also computes the Fundamentals using a harmonic pertebation, which can be compared to the actual analysis.
In this way you can verify that all of the transformations are being computed correctly. 
This code only works for the last 3 dimensions defined for Vmax (hard-coded the potential analysis). 

# Associated Files

## Compile.sh
A simple bash script for running the code. Simply run the following command on terminal to execture; `$ bash compile.sh`.

## qmc_pot_test.f90
The program main.

## Input.dat
A sample input file for defining the excitations, input geometry, number of evaluations and etc.
Check the main program to see the actual read statements.

## ssobol_stdnormal.f90
Fortran program to transform sobol points (0,1) uniformly distributed, to a normal distribution.
We utilize the Beasley-Springer-Moro algorithm to complete this transformation. 

## Running The Code
To run this code you will need an input geometry (xyz) for a water monomer. 
The code has been simplified to run using the TIP4P potential surface (only), and needs the TIP4P Fortran Module to run.
The code can only analyze a system with the last 3 dimesnions containing pertebations (Vmax = 0,0,0,0,0,0,#,#,#)
It also required the Fortran scrambled sobol point generator provided by the MCQMC wiki page. 
Finally we use the cg.f Fortran module for all matrix operations. 

## Analysis Example:
To verify the method, I provide an example calculation.
For a pertebation of 0.02 take the value from the fundamentals.dat file (fund) and compare to the analagous value in the harmonic.dat file (harmonic).

fund - ((fund - (fund / 1.02))/2) = harmonic
