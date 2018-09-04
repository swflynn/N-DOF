# 1D_HO_matrix
Fortran 90 Implementation for evaluating the Potential Energy Matrix for a Quantum Harmonic Oscillator (Hermite Polynomials) for a single spatial dimension. 
The Potential Energy Matrix elements are computed numerically using Quasi-Monte Carlo methods (both standard and scrambled sobol sequences).  
The Basis Set is orthonormal, therefore off-diagonal = 0,diagonal = 1. 
The user provides what excitation to compute up to, and the number of sobol points to use for QMC.
Matrix elements are composed of all combinatorics in the excitations. 
Sobol points are generated over a uniform distribution, and converted to a normal distribution using the Beasley-Springer-Moro algorithm for QMC. 

## Standard
Uses the sobol.f90 module to generate a standard sobol sequence for QMC (see FSU John Burkardt for sobol.f90 source code). 

sobol_stdnormal.f90 ==> 

Transforms sobol points (0,1) uniformly distributed through the Beasley-Springer-Moro algorithm to a normal distribution. 

## Scrambled
Requires an input sequence to be used for QMC (assumes input sequence is a scrambled sobol sequence). 
An example script for generating scrambled sobol sequences using Matlab is provided. 

s_sobol.m ==>

Matlab script for generating scrambled sobol sequences and storing to file 

ss_unif ==>

An example scrambled sobol sequence stroed to file

ss_stdnormal.f90 ==> 

Takes in a vector of sobol points (0,1) uniformly distributed, and uses the Beasley-Springer-Moro algorithm to transform them to a normal distribution. 

## Data_Analysis
Python 2 script for plotting various potential the convergence as a function of iteration for the Potential Energy Matrix Elements. 
