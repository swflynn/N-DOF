% Script to generate scrambled sobol sequence in matlab
%Determine how long it takes to generate points
ti = cputime;
%define points, skip, dimensionality of the sequence
Nsobol = 1;
Nsobol_skip = 1;
d = 2;

% generate the sequence and save it to a file
points = (sobolset(d, 'Skip', Nsobol_skip));
points = (scramble(points, 'MatousekAffineOwen'));
sequence = double(net(points, Nsobol));
dlmwrite('s_sobol_unif.dat', sequence , 'delimiter', ' ', 'precision', 8);
% write to terminal how many points and how long it took
tf = cputime;
sprintf('The total time to write out all the points was')
t = tf - ti
sprintf('We had this many sobol points')
Nsobol
sprintf('We had this many spatial dimensions')
d
