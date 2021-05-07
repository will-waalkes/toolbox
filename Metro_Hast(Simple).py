# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:49:30 2017

@author: willwaalkes
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt

def f_x(X):
    return 1./(1.+X**4.) #Our model PDF - this turns into the Posterior probability

N = 100000 #Number of steps
x_0 = 500 #Initial Guess

x = np.zeros(shape=(N+1)) #Empty arrays to fill with updated x or f(x) values
fx = np.zeros(shape=(N+1))
x[0] = x_0
fx[0] = f_x(x_0)
for i in range(N):
    x_i = np.random.normal(x_0,3) #Choose a random number from a Gaussian
    u_i = np.random.uniform(0,1) #Choose an acceptance value

    alpha = np.amin([1,(f_x(x_i)/f_x(x_0))],axis = 0)

    if u_i < alpha:
        x[i+1] = x_i
        fx[i+1] = f_x(x_i)
        x_0 = x_i
    else:
        x[i+1] = x_0
        fx[i+1] = f_x(x_0)
        x_0 = x_0
    
plt.plot(x)
plt.ylabel('X')
plt.xlabel('N')
plt.title('Metropolis-Hastings')
plt.show()
plt.clf()

#Make a histogram 
plt.hist(x[50000:-1], 100, normed=True)

x_vals = np.linspace(-10,10,81)
func = 0.45*f_x(x_vals)
plt.plot(x_vals,func,color='r')
plt.xlabel('Sampled X')
plt.ylabel('Normalized Probability Density')
plt.title('Metropolis Hastings Reults with PDF')
plt.yscale('log')
