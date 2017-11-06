10/25/17
Corrected version of the code. 
We found an error in the transformations, we were scaling the harmonic Hamiltonian into normal coordinatesd with an omega^2 it should be onega
(note hbar defined the energy, the omega cancels to make it unitless). 

This version of the code corrects for this problem. 
We are also missing a factor of 2 in the code, that we need to find

We are more confident in the code now, because we compute the harmonic potential using both the transformaiton into normal modes, and using the tradsitional definition with the Hessian. 
The only difference between the results si teh factor of 2 we are missing. 

This means we are (almost) computing the correct potnetial with our transformations. 
