import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('data_file.txt')

nsobol = data[:,][:,0]

diag_elem = data[:,][:,1:5]

diag_elem = abs(1-diag_elem)

#print diag_elem[0:5]

off_diag_1 = data[:,][:,5:9]
off_diag_1 = abs(off_diag_1)

off_diag_2 = data[:,][:,9:13]
off_diag_2 = abs(off_diag_2)

#plot off diagonal 
plt.figure(1)
plt.plot(nsobol,diag_elem)
plt.legend(['<3|3>', '<5|5>', '<7|7>', '<9|9>'])
plt.xlabel('Number of Sobol Points')
plt.ylabel('Relative Error')
plt.savefig('Diag_Conv_10E8.png')
#plt.show()


plt.figure(2)
plt.plot(nsobol,off_diag_1)
plt.legend(['<0|5>', '<0|9>', '<5|0>', '<9|0>'])
plt.xlabel('Number of Sobol Points')
plt.ylabel('Distance From 0')
plt.savefig('Off_Diag_Conv_10E8_1.png')
#plt.show()

plt.figure(3)
plt.plot(nsobol,off_diag_1)
plt.legend(['<6|2>', '<5|7>', '<7|9>', '<8|3>'])
plt.xlabel('Number of Sobol Points')
plt.ylabel('Distance From 0')
plt.savefig('Off_Diag_Conv_10E8_2.png')
#plt.show()

