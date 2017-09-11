# Vibrational Spectroscopy; Quasi-MC Methods
Fortran-90 implementation of N-dimensional (N<10), k excitation (k<10) vibrational spectroscopy analysis for water clusters.

Quasi Monte Carlo methods for calculating a potential energy difference matrix to analyze anharmonic corrections associated with the standard harmonic approximation.
Utlilizes both standard sobol and scrambled sobol quasi-random sequences for evaluating high-dimensional integrals in the potential energy matrix.
Matrix Element Convergence and Eigenvalue Convergence as a function of point available. 

### Current Developments
Adding a potential energy surface for water (TIP4P and MB-POL), and utilizing a normal mode basis for analyizing a single water molecule.
Future work will apply this methodology combined with a local monomer approximation (a single water molecule as the unit) to study water clusters. 

## Author
Shane W. Flynn, Vladimir A. Mandelshtam. 2017. UCI

## References
#### General Method and Application Overview
http://aip.scitation.org/doi/abs/10.1063/1.4829836 

#### Numerical Methods
http://aip.scitation.org/doi/abs/10.1063/1.4788977

http://aip.scitation.org/doi/abs/10.1063/1.4754819

#### Local Monomer Approximation
http://pubs.acs.org/doi/abs/10.1021/jp5061182

#### Potential Energy Surfaces
http://paesanigroup.ucsd.edu/index.html

http://pubs.acs.org/doi/pdf/10.1021/ct400863t

#### Quasi MC Code:
Currently utilizing the scrambled sobol Fortran code made available on the MCQMC wiki page.
http://roth.cs.kuleuven.be/wiki/Public_software

Regards to John Burkardt (FSU) for his advice on scrambling algorithms and Fortran resources. 
