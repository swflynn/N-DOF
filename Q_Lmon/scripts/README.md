# Scripts
Collection of various programs written during project development.

### 1D_HO_matrix:
Potential Energy Matrix evaluated numerically (Quasi-Monte Carlo) for a single dimension Quatum Harmonic Oscillator.
User defines highest degree Hermite Polynomial (matrix is composed of all combinatorics over these excitations), and uses Normally Distributed Sobol Points for quasi-Monte Carlo integration. 

### nD_PE_matrix:
Fortran-90 script for calculating any dimensional (D<10) any excitation (k<10) Harmonic Oscillator Potential Energy matrix.
Modified hermite polynomials are used to compute the wavefunction, and QMC Sobol Sequences (normally distributed, Beasley-Springer-Moro) are used to perform numerical integration.
Eigenvalue and matrix element convergance as a function of iteration are available. 
No external field exists, wavefunctions are orthonormal, therefore matrix elements should be either 1(diagonal), or 0, eigenvalues should be 1. 
Sobol sequences are generated using `Sobol.f90`.

### nD_s_PE_matrix:
Fortran-90 script for calculating any dimensional (D<10) any excitation (k<10) Harmonic Oscillator Potential Energy matrix.
Modified hermite polynomials are used to compute the wavefunction, and QMC Scrambled-Sobol Sequences (normally distributed, Beasley-Springer-Moro) are used to perform numerical integration.
Eigenvalue and matrix element convergance as a function of iteration are available. 
No external field exists, wavefunctions are orthonormal, therefore matrix elements should be either 1(diagonal), or 0, eigenvalues should be 1. 
Scrambled Sobol Sequences are generated using the matlab `scramble` utility, part of the `qrandset` class.
The code simply needs a datafile containing each point to be evaluated.

### 1D_s_HO_matrix:
Fortran-90 script for calculating a 1D Harmonic Oscilator matrix (weighted by a Gaussian).
The code allows users to define which degree polynomial to calculate up to for the HO Matrix elements, and uses Normally Distributed Scrambled Sobol Points for quasi-Monte Carlo integration. 
The points must be generated externally (a matlab script has been attached for doing this) see
https://www.mathworks.com/help/stats/qrandset.scramble.html for matlab details. 

### permutation:
Fortran-90 script to determine the number of permutations available in a system with a d-dimensional wavefunction and a maximum excitation constraing of Vmax (the total number of excitations allowed in the system).

### PES_test:
Fortran-90 script that evaluates the potential energy directly from the Hessian, and using the normal mode transformations.
This test code is stripped down, and test to be sure all of hte normal mode transformation are being done accurately.
Each method should produce the same result.

### Harmonic:
Fortran-90 script that computes the Fundamental Frequencies for a water cluster using the TIP4P potential only.
This code is a stripped-down version that also computes the fundamental frequencies using a harmonic pertebation.
In this way you can determine directly if all transformations are appropriately defined.

### nD_general_PEmat
Generalized code using ssobol.f 
Code only computes the PE matrix was the previous main need to make this github write-up

### qMC_nm:
Need to make this github ready. First code for monomer doing fundamental analysis. 
(Currently the program mian as of 4-27-18). 

# Open Source Code
The following codes are/have been used during the project.
Please see the original documentation for more details, and give proper reference to the authors. 
Many thanks to all the authors who have made these codes available, I have found them to be extremly helpful throughout this project.

## Sobol.f90 
The Fortran-90 module used for Sobol Sequence Generation.
Taken from http://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html

## ACM Collected Algorithms
A collection of various algorithms, notably the various quasi-random sequence generators (647, 659)
Please see http://calgo.acm.org/ for source code. 

## MCQMC Wiki Page
Public Software containing MC, QMC, MCMC software for the community. 
Notably the `ssobol.f` fortran source code has been utilized for scrambled sobol sequences for MC integration. 
Please see http://roth.cs.kuleuven.be/wiki/Main_Page for the wiki page, and various resources available.
