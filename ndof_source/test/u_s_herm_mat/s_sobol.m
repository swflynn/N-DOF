% Script to generate scrambled sobol sequence in matlab

%define points, skip, dimensionality of the sequence
Nsobol = 100000;
Nsobol_skip = 10000;
d = 1;

% generate the sequence and save it to a file
points = (sobolset(d, 'Skip', Nsobol_skip));
points = (scramble(points, 'MatousekAffineOwen'));
sequence = double(net(points, Nsobol));
dlmwrite('s_sobol_unif.dat', sequence , 'delimiter', ' ', 'precision', 15);
