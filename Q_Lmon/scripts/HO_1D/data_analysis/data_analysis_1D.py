#=============================================================================80
#                      	    HO 1D Data Analysis
#=============================================================================80
#    		Discussion:
# Python 2 implementation for plotting convergence as a function of iteration 
# for the 1D Harmonic Oscillator Quasi-Monte Carlo Potential Energy Matix.
# Plot convergence and absolute error for diagonal elements.
# Plot convergence as a function of iteration for all matrix elements.
#    		Modified:
# 20 March 2017
#    		Author:
# Shane Flynn
#==============================================================================#
from sys import argv
import pandas as pd                     
import numpy as np                     
import matplotlib.pyplot as plt       
#==============================================================================#
#                               Variables
# script        ==> data_analysis_1D.py
# data_file     ==> Convergence as a function of iteration
# exc           ==> Maximum excitation used in simultaion
# x             ==> Iteration
# diag_index    ==> Diagonal Matrix Elements (Index)
#==============================================================================#
#			 Read Simulation Parameters
# See run.sh for syntax
#==============================================================================#
script, data_file, exc = argv
exc = int(exc)
#==============================================================================#
#			    Read Convergence Data
# Iteration is the first element of every row
#==============================================================================#
df_conv = pd.read_csv(data_file, delim_whitespace=True,header=None)
x = df_conv[0]
#==============================================================================#
#                 Compute location of diagonal matrix elements
#==============================================================================#
diag_index = []
index=0
counter=0
while index < exc**2:
    index = counter*exc + counter + 1
    diag_index.append(index)
    counter += 1
print 'Index of Diagonal Elements'
print diag_index
#==============================================================================#
#                     Plot Diagonal Elements Convergence
#==============================================================================#
for i in diag_index:
    plt.plot(x,df_conv[i], label='%s'%i) 
plt.xlabel('Nsobol')
plt.ylabel('Convergence (Diagonal)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) # legend outside 
plt.savefig('diag_conv.png', bbox_inches='tight') 
plt.close()
#==============================================================================#
#                       Plot Absolute Error Diagonal
#==============================================================================#
for i in diag_index:
    plt.plot(x,abs(1-df_conv[i]), label='%s'%i)
plt.xlabel('Nsobol')
plt.ylabel('Absolute Error (Diagonal)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
plt.savefig('diag_conv_abs.png', bbox_inches='tight')
plt.close()
#==============================================================================#
#                   Plot Absolute Error (Yaxis = Log Scale) 
#==============================================================================#
for i in diag_index:
    plt.plot(x,abs(1-df_conv[i]), label='%s'%i)
plt.xlabel('Nsobol')
plt.ylabel('Absolute Error (Log Scale)')
plt.yscale('log')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
plt.savefig('diag_log_conv_abs.png', bbox_inches='tight')
plt.close()
#==============================================================================#
#                      All Matrix Elements Convergence 
#==============================================================================#
i = 1
while i < exc**2+1:
    plt.plot(x,df_conv[i])
    plt.xlabel('Nsobol')
    plt.ylabel('Potential Energy')
    plt.savefig('1D_conv_all_%s.png'%i) 
    plt.close()
    i += 1
