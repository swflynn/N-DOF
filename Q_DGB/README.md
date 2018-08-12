# Q_DGB
QuasiRandom Distrbiuted Gaussian Basis with QuasiMonteCarlo

## Overview
Fortran-90 implementation of Quasirandom Gaussian Basis Sets combined with Quasi Random Monte Carlo for evaluating water Clusters. 

The program applies Quasi Monte Carlo methods for numerical integration of anharmonic terms in the fundamental frequencies and subsequent vibrational spectra of a water cluster. 
Currently Sobol and Scrambled Sobol Sequences can be used to compute the normal modes, and fundamental frequencies for a sigle water atom within a water cluster. 

### Current Developments
Make a draft of the code, generte Multi-Variant Gaussian Locations and Compute Potential Matricies at each Gaussian. 
Start with a simple Potential that we can solve analytically (multi-dimensional quadratic potential). 

## References
A small collection of useful papers and content motivating the project. 

#### General Method and Application Overview
Gaussian Bases and their sampling distributions are well established in the literature, for example [Garashchuk and Light](https://aip.scitation.org/doi/abs/10.1063/1.1348022). 

The application of quasi-random sequences to compute ground state energies of quantum systems can be found in the literature; for example [Georgescu and Mandelshtam](https://aip.scitation.org/doi/abs/10.1063/1.4829836).

#### Quasirandom Sequences:
Currently utilizing Fortran Code for the generation of Sobol Sequences made available by [John Burkardt et. al](https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html).

The Monte Carlo and Quasi-Monte Carlo Wiki page ([MCQMC](http://roth.cs.kuleuven.be/wiki/Main_Page)) is another useful resource for generating and applying Monte Carlo methods. 

## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2018. UCI
