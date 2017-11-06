11/1/17
This version of the code correctly calls the potential
We verify by calling the potential with our normal mode transformations, and with the hessian
We have corrected the factor of 2 issue that was present before!

We now need to do two things
1. Replace the difference potential with a hamronic potential*0.1
this is just as mall perterbation of the harmonic potential, the results should be extremlyclose to the harmonic potential which are well known. 

2. Once we get that code working we can compute the true difference potential and determine the fundamental frequencies, and hopefully everything will work out!
