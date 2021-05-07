# -*- coding: utf-8 -*-
"""
Created on Sat May 20 12:11:58 2017

@author: willwaalkes
"""

from random import random,randrange
from math import exp,pi
from numpy import ones
import os
import sys
import batman
import numpy as np
import matplotlib.pyplot as plt

T = 10.0 #Temperature
N = 1000 #number of particles
steps = 250000

#create a 2D array to store the quantum numbers
n = ones([N,3],int)

#Main Loop
eplot = []
E = 3*N*pi*pi/2 #Simplified energy equation - this is the model.

for k in range(steps):
    #choose the particle and the move
    i = randrange(N)
    j = randrange(3)
    if random()<0.5:
        dn = 1
        dE = (2*n[i,j]+1)*pi*pi/2
    else:
        dn = -1
        dE = (-2*n[i,j]+1)*pi*pi/2
        
    #decide whether to accept the move
    if n[i,j]>1 or dn==1:
        if random()<exp(-dE/T):
            n[i,j] += dn
            E += dE
        
    eplot.append(E)

# Make the graph
plt.plot(eplot)
plt.ylabel("Energy")
plt.show()

os.sys.exit()
# NOW WITH BATMAN

#Next we create a TransitParams object to store the physical parameters describing the transit:

params = batman.TransitParams()
params.t0 = 0.                       #time of inferior conjunction
params.per = 1.628930                #orbital period
params.rp = 0.0512                   #planet radius (in units of stellar radii)
params.a = 16.                       #semi-major axis (in units of stellar radii)
params.inc = 90.                     #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1956, 0.3700]          #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

#Note that for circular orbits, batman uses the convention params.w = 90. 
#The units for params.t0 and params.per can be anything as long as they are consistent.

#We also need to specify the times at which we wish to calculate the model:

t = np.linspace(-0.05, 0.05, 500)

#Using these parameters, we initialize the model and calculate a model light curve:

m = batman.TransitModel(params, t)    #initializes model

for i in range(len(N)):
    flux = m.light_curve(params)          #calculates light curve

    #Now simulate data points

    sim = flux + np.random.normal(0,0.02,len(t))*np.sqrt(N[i])#/N[i]

    plt.plot(t, flux)
    plt.scatter(t, sim)
    plt.xlabel("Time from central transit")
    plt.ylabel("Relative flux")
    plt.title('Counts per minute (N):')
    plt.title(N[i], loc = 'right')
    plt.show()
    plt.clf()
'''
steps = 250000

#Main Loop
eplot = []

for k in range(steps):
    #choose the particle and the move
    if random()<0.5:

    else:
        
    #decide whether to accept the move
    if random()<: # NEED A PROBABILITY CALCULATION
        R += dR
        
    eplot.append(R)

# Make the graph
plt.plot(eplot)
plt.ylabel("Planet Radius")
plt.show()
'''