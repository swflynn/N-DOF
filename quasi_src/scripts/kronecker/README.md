# Hamiltonian Calculation
H = \delta_vv * \epsilon_v + Potential Difference Matrix

This code calculates the kronecker delta times the energy, both of which are a function of the v, the permutations. 
I will refer to this matrix as `kd_mat`. 

\epsilon_v = \sum_k=1^d hbar omega_k (\nu_k + 0.5)

## Matrix details
We are adding this matrix to the potential difference, therefore it is size  kd_mat(Jmax,Jmax)

Jmax is the total number of permutations available given a spatial dimension and a maximum excitation (Vmax).

For the example code: spatial_dim = 3, Vmax = 2, there are therefore 10 permutations in total. 
000, 001, 002, 010, 011, 020, 100, 101, 110, 200 

The quantum number is given by the sum of each excitation in each spatial dimension therefore 

0 , 1 , 2, 1, 2, 2, 1, 2, 2, 2

## More on the Matrix
The matrix is Jmax by Jmax, each permutation set defines the rows and columns.
Because the kd, off-diagonal elements must be 0, therefore the diagonal elements must be calcualted. 
These diagonal elements depend on the quantum number, each diagonal has a set of excitations, adding each spatial dimension contribution gives the quantum number. 

## The Code
The code takes in a vector of omega values (the eigenvalues, size = spatial dimension). 
It then takes in a given Vmax, Spatial_dimension to give Jmax, and v(d,Jmax) (the set of permutations). 

We then set kd_mat(Jmax,Jmax) the same dimension as the potential energy difference matrix. 
The off-diagonal elements are set to 0 due to the kd, the permutations are different in each spatial dimension, therefore the product is 0. 
For example take the first row, second column element. 
000, 001. If you evaluate these, the first to spatial dimensions would give a value of 1, but the last dimension is a different function and gives a value of 0, therefore the product is 0. 

We then populate our diagonal matrix elements, by determining the quantum number of that specific permutation. 
We then use the equation to get the energy associated by summing all terms up to the excitation.
There may be a mismatch here, quantum number starts at 0, but omega starts at 1, so this may be a potential issue to check out. 
