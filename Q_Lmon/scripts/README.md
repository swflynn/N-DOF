# Scripts
Collection of various programs written during project development.

### HO 1D:
Potential Energy Matrix evaluated numerically (Quasi-Monte Carlo) for a single dimension Quatum Harmonic Oscillator.
User defines highest degree Hermite Polynomial (matrix is composed of all combinatorics over these excitations), and uses Normally Distributed Sobol Points for quasi-Monte Carlo integration. 

### HO nD:
(Standard)
Fortran-90 script for calculating any dimensional (D<10) any excitation (k<10) Harmonic Oscillator Potential Energy matrix.
Modified hermite polynomials are used to compute the wavefunction, and QMC Sobol Sequences (normally distributed, Beasley-Springer-Moro) are used to perform numerical integration.
Eigenvalue and matrix element convergance as a function of iteration are available. 
No external field exists, wavefunctions are orthonormal, therefore matrix elements should be either 1(diagonal), or 0, eigenvalues should be 1. 
Sobol sequences are generated using `Sobol.f90`.

(Scrambled)
Fortran-90 script for calculating any dimensional (D<10) any excitation (k<10) Harmonic Oscillator Potential Energy matrix.
Modified hermite polynomials are used to compute the wavefunction, and QMC Scrambled-Sobol Sequences (normally distributed, Beasley-Springer-Moro) are used to perform numerical integration.
Eigenvalue and matrix element convergance as a function of iteration are available. 
No external field exists, wavefunctions are orthonormal, therefore matrix elements should be either 1(diagonal), or 0, eigenvalues should be 1. 
Scrambled Sobol Sequences are generated using the matlab `scramble` utility, part of the `qrandset` class.
The code simply needs a datafile containing each point to be evaluated.

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

# Open Source Codes
The following codes have been useful during this project.
Refer to the original source code and documentation for more details.
Many thanks to all the authors who have made these codes available (special thanks to John Burkardt for his knowledge on quasi-random sequences and their applications). 

## Sobol.f90 
Fortran-90 module for generating standard Sobol Sequences. 

http://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html

## ACM Collected Algorithms
A collection of various algorithms, notably the various quasi-random sequence generators (647, 659)

http://calgo.acm.org/ for source code. 

## MCQMC Wiki Page
Public Software containing MC, QMC, MCMC programs.

http://roth.cs.kuleuven.be/wiki/Main_Page 
