# Vibrational Spectroscopy; Quasi-MC Methods
Fortran-90 implementation of N-dimensional (N<10), k excitation (k<10) vibrational spectroscopy analysis for water clusters.

The program applies Quasi Monte Carlo methods for numerical integration of anharmonic terms in the fundamental frequencies and subsequent vibrational spectra of a water cluster. 
Currently Sobol and Scrambled Sobol Sequences can be used to compute the normal modes, and fundamental frequencies for a sigle water atom within a water cluster. 

### Current Developments
Compute the fundamental frequencies for every water monomer in the cluster as indepedent calculaitons.
Automate the data analysis for every calculation. 

## References
The motivation for this project comes in part from previous works suh as the Multimode Software by Bowman and etc. 
Some reference papers for context are:

#### General Method and Application Overview
http://aip.scitation.org/doi/abs/10.1063/1.4829836 

#### Numerical Methods
http://aip.scitation.org/doi/abs/10.1063/1.4788977

http://aip.scitation.org/doi/abs/10.1063/1.4754819

#### Local Monomer Approximation
http://pubs.acs.org/doi/abs/10.1021/jp5061182

#### Potential Energy Surfaces
The MBPOL PES is currently being used for these calculations.
http://paesanigroup.ucsd.edu/index.html
http://pubs.acs.org/doi/pdf/10.1021/ct400863t

TIP4P PES has also been implemented. 

#### Quasi MC Code:
Currently utilizing the scrambled sobol Fortran code made available on the MCQMC wiki page (ssobol.f).
http://roth.cs.kuleuven.be/wiki/Public_software

## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2017. UCI
