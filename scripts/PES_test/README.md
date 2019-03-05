# PES_test
Code for evaluating the Potential Energy Surface using two different methods.
It can be calculated directly from the Hessian, or through normal mode analysis.
The project utilizes normal modes to make the Hamiltonian Seperable, and this script provides a test to ensure our transfomrations are correctly defined (the two calculations should provide the exact same answer).

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
It also required the Fortran scrambled sobol point generator provided by the MCQMC wiki page. 
Finally we use the cg.f Fortran module for all matrix operations. 
