#!/usr/bin/env python
#=============================================================================80
#                       Gaussian Distributed bins
#=============================================================================80
#       Discussion: 
# Python 2 implementation for computing 1D gaussian bins
#==============================================================================#
#       Modified:
# 13 September 2018
#       Author:
# Shane Flynn
#==============================================================================#
import numpy as np
from sys import argv
from scipy.special import erfinv
from scipy.special import erf
#==============================================================================#
#==============================================================================#
# argv = script, lower_bound, upper_bound, Number_Bins
#==============================================================================#
script,a,b,Npoints = argv
a = float(a)
b = float(b)
Npoints = int(Npoints)
#==============================================================================#
# Compute the area from b,a under the gaussain
# int dx e^{-x^2/2} = sqrt(pi/2)erf(x/sqrt(2))
#==============================================================================#
area = np.sqrt(np.pi/2)*(erf(b/np.sqrt(2.))-erf(a/np.sqrt(2.)))
#print('Area under Gaussian [a,b] ==> %r')% area
#==============================================================================#
# Compute width to define each bin
#==============================================================================#
bins = area / Npoints
#print('Area/ Bin ==> %r')% bins
#==============================================================================#
# define starting point x_0 = 0
#==============================================================================#
x0 = 0
x = []
x.append(x0)
#print x
#print type(x[0])
#x1 = np.sqrt(2)* erfinv( bins / (np.sqrt(np.pi/2)) + erf(x0 / np.sqrt(2))  )
#print('x1 ==> %r')% x1
#==============================================================================#
# Iterate over [a,b]
# Python does not incluse upper bound, so add 1
#==============================================================================#
for i in range(1,(Npoints+1)/2):
    xi = np.sqrt(2)* erfinv( bins / (np.sqrt(np.pi/2)) + erf(x[i-1] / np.sqrt(2))  )
    x.append(xi)
#print x
#==============================================================================#
# Gaussian is symetic, choose negatice points
#==============================================================================#
y = []
y.append(-x0)
for i in range(1,(Npoints+1)/2):
    yi = -x[i]
    y.append(yi)
#print y
#==============================================================================#
# Remove initial point from y to avoid redundancy
#==============================================================================#
y.remove(y[0])
#print y
#==============================================================================#
# merge lists for convenience
#==============================================================================#
data = x + y
#print data
#==============================================================================#
# Compute all permutations and write to file
#==============================================================================#
k=0
for i in range(0,(Npoints)):
    for j in range(0,(Npoints)):
        k +=1
        print('%r   %r')%(data[i],data[j])
#print k
#print Npoints*Npoints
