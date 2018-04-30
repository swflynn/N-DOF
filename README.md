# Vibrational Spectroscopy; Quasi-MC Methods
Fortran-90 implementation of N-dimensional (N<10), k excitation (k<10) Fundamental Frequency analysis for water clusters.

The program applies Quasi Monte Carlo methods for a numerical framework to compute the RoVibrational Hamiltonian of small water clusters.
Ultimately Fundamental Frequencies for a given water cluster and Potential Energy Surface are Computed using normal mode coordinates and the local monomer approximation, combined with Sobol or Scrambled Sobol Quasi Monte Carlo methods. 

## References
RoVibrational Energy Calculations are an essential goal of Computational Chemistry. 
Some related references and useful software are provided below for context. 

#### General Method and Application Overview
Quasi-Random Monte Carlo Methods applied to small water clusters can be found in the literature, for example see the work of [Georgescu and Mandelshtam](http://aip.scitation.org/doi/abs/10.1063/1.4829836).

Basis construction and numemrical methods for computing fundamental frequencies through QMC can be found in the literature, for example see [Brown and Mandelshtam](http://aip.scitation.org/doi/abs/10.1063/1.4788977).

The application of the Local Monomer Approximation in the context of RoVibrational Spectroscopy can be found uin the Literature, see for example [Liu and Bowman](http://pubs.acs.org/doi/abs/10.1021/jp5061182).

#### Quasirandom Sequences:
The Monte Carlo and Quasi-Monte Carlo Wiki page ([MCQMC](http://roth.cs.kuleuven.be/wiki/Main_Page)) is a useful resource for generating and applying Monte Carlo methods. 
We are currently utilizing the [Scrambled Sobol Sequence Generator](http://roth.cs.kuleuven.be/wiki/Public_software) made available for Fortran implementations (ssobol.f). 

## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2017. UCI

