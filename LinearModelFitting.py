# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:14:06 2017

@author: willwaalkes
"""

import os
import sys
import batman
import numpy as np
from numpy.linalg import inv
from numpy.linalg import eig
from numpy import dot
import matplotlib.pyplot as plt

# Make a basic y = mx + b fit from David Hogg's linear model data.

'''Y = np.matrix([[592.0],[401.0],[583.0],[402.0],[495.0],[173.0],[479.0],[504.0],
               [510.0],[416.0],[393.0],[442.0],[317.0],
               [311.0],[400.0],[337.0],[423.0],[334.0],[533.0],[344.0]])
A = np.matrix([[1.,201.],[1.,244.],[1.,47.],[1.,287.],[1.,203.],[1.,58.],[1.,210.],
               [1.,202.],[1.,198.],[1.,158.],[1.,165.],
               [1.,201.],[1.,157.],[1.,131.],[1.,166.],[1.,160.],[1.,186.],[1.,125.],
                [1.,218.],[1.,146.]])
A_Plot = [201,244,47,287,203,58,210,202,198,158,165,201,157,131,166,160,186,125,218,146]
Y_Plot = [592,401,583,402,495.0,173.0,479.0,504.0,510.0,416.0,393.0,442.0,317.0,
               311.0,400.0,337.0,423.0,334.0,533.0,344.0]
'''               
# PART 3  Add another column to matrix A containing the values x^2 and another element to vector X (call it q).
Y = np.matrix([[495.0],[173.0],[479.0],[504.0],
               [510.0],[416.0],[393.0],[442.0],[317.0],
               [311.0],[400.0],[337.0],[423.0],[334.0],[533.0],[344.0]])
A = np.matrix([[1.,203.,203.**2],[1.,58.,58.**2],[1.,210.,201.**2],
               [1.,202.,202.**2],[1.,198.,198.**2],[1.,158.,158.**2],[1.,165.,165.**2],
               [1.,201.,201.**2],[1.,157.,157.**2],[1.,131.,131.**2],[1.,166.,166.**2],
                [1.,160.,160.**2],[1.,186.,186.**2],[1.,125.,125.**2],
                [1.,218.],[1.,146.]])
A_Plot = [203,58,210,202,198,158,165,201,157,131,166,160,186,125,218,146]
Y_Plot = [495.0,173.0,479.0,504.0,510.0,416.0,393.0,442.0,317.0,
               311.0,400.0,337.0,423.0,334.0,533.0,344.0]

               
sigy = np.array([21.,15.,27.,14.,30.,16.,14.,25.,52.,16.,34.,31.,42.,26.,16.,22.])
sigy = sigy**2
sigx = [5,9,4,4,11,7,5,5,5,6,6,5,9,8,6,5]
k = len(sigy)
C = np.zeros((k,k,k))


for i in range(k):
    for j in range(k):
        for l in range(k):
            if i == j == l:
                C[i][j][l] = sigy[i]
            else: C[i][j][l] = 0
        

cov = inv([A.transpose()*inv(C)*A])
z = [A.transpose()*inv(C)*Y]
X = dot(cov,z)

print cov
print X

os.sys.exit()

b = float(X[0][0])
m = float(X[0][1])
q = float(X[0][2])

x = np.linspace(0,300,100)
y = m*x+b

#print A, A_Plot
#os.sys.exit()

plt.plot(x,y)
plt.scatter(A_Plot,Y_Plot)
plt.errorbar(A_Plot,Y_Plot, yerr=np.sqrt(sigy), xerr=None, linestyle ="None")
plt.ylim(0,700)
plt.xlim(0,300)
plt.title("Unnecessary and wrong linear fit")
plt.show()