# Permutations With Parameters
As with the other permutations script we are considering a systems with a few paramaters. 
There are the number of basis functions you are going to use (the number of spatial dimensions), the excitation associted with each basis function, and teh total excitation of the system. 

The combinatorics scaling of the problem makes this scale very rapidly, so we have developed a different permutation algorithm allowing for more control. 

# Parameters
To use this code you need to define the following:

Dim :: Integer; this sets the number of spatial dimensions available in the system (9 is the maximum value). 
To run this code you must keep Dim=9 as we are studying a single water molecule, and our ?Hessian will be 9 by 9. 

Vtot :: Integer; This sets the total excitation allowed for the entire system (meaning the sum of the excitations across each spatial dimension). 

Vmax(9) :: Integer; This sets the maximum excitation allowed in any one spatial dimension. 
We have hard coded this system to study a water monomer, therefore Vmax(9) is of size 9 for the 3N dimensions associated with a single water molecule. 

Vtot restricts the combinations available for given Vmax. 

You MUST provide a Vmax vector for the code i.e. Vmax = [0,0,0,0,0,0,3,3,3]

v(Jmax,Jmax) :: Integer, allocatable; contians all of the combinations available given Vmax, Vtot, Dim. 
