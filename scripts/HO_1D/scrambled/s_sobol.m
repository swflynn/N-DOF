%=============================================================================80
%		      Scrambled Sobol Sequence (Matlab)
% Matlab Quasi-Random sequence generator
%=============================================================================80
%           Discussion:
% Matlab Quasi-Random sequence generator, written to file. 
% Uses the MatousekAffineOwen scrambling method for sobol sequences.
%==============================================================================%
% See Matlab Documentation:
% https://www.mathworks.com/help/stats/sobolset.html
%==============================================================================%
%           Modified:
% 20 March 2017
%           Author:
% Shane Flynn
%==============================================================================%
%                               Variables
%==============================================================================%
% Nsobol        ==> Number of Scrambled Sobol Points
% Nsobol_skip   ==> Skip for Generator (suggested=Nsobol)
% dimen         ==> Dimensionality of Sequence
%==============================================================================%
%                           Define Sequence 
%==============================================================================%
Nsobol = 100;
Nsobol_skip = 100;
dimen = 1;
%==============================================================================%
%                       Generate and Write to File
%==============================================================================%
points = (sobolset(dimen, 'Skip', Nsobol_skip));
points = (scramble(points, 'MatousekAffineOwen'));
sequence = double(net(points, Nsobol));
dlmwrite('s_sobol_unif.dat', sequence , 'delimiter', ' ', 'precision', 15);
