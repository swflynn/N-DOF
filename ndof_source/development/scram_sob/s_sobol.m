% Script to generate scrambled sobol sequence in matlab

%define points, skip, dimensionality of the sequence
t_ini = cputime
Nsobol = 100000;
Nsobol_skip = 100000;
d = 3;

% generate the sequence and save it to a file
points = (sobolset(d, 'Skip', Nsobol_skip));
points = (scramble(points, 'MatousekAffineOwen'));
sequence = double(net(points, Nsobol));
dlmwrite('s_sobol_unif.dat', sequence , 'delimiter', ' ', 'precision', 15);
t_fin = cputime
time = t_fin - t_ini