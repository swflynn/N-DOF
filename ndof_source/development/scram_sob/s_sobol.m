% Script to generate scrambled sobol sequence in matlab

%define points, skip, dimensionality of the sequence
Nsobol = 10;
Nsobol_skip = 10;
d = 1;

% generate the sequence and save it to a file
points = sobolset(d, 'Skip', Nsobol_skip);
points = scramble(points, 'MatousekAffineOwen');
sequence = net(points, Nsobol); 
dlmwrite('s_sobol.dat', sequence , ' ');
