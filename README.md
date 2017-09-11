# Vibrational Spectroscopy; Quasi-MC Methods
Fortran-90 implementation of N-dimensional (N<10), k excitation (k<10) vibrational spectroscopy analysis for water clusters.

The code utilizes quasi Monte Carlo methods for calculating a potential enegy difference matrix to analyze anharmonic corrections associated with the standard harmonic approximation.
The code currently utlilizes both standard sobol and scrambled sobol sequences for evaluating high-dimensional integrals in the potential energy matrix.
Matrix Elements Convergence and Eigenvalue Convergence available. 

## Current Developments
Currently working on adding a potential energy surface for water, and utilizing a normal mode basis for analyizing a single water molecule.

Future work will apply this methodology combined with a local monomer approximation (a single water molecule as the unit) to study water clusters. 

## Author:
Shane Flynn, Vladimir Mandelshtam. 2017. UCI

### References:
For motivation of the method see:
http://aip.scitation.org/doi/abs/10.1063/1.4829836 

For motivation of the basis and numerical techniques see:
http://aip.scitation.org/doi/abs/10.1063/1.4788977
http://aip.scitation.org/doi/abs/10.1063/1.4754819
