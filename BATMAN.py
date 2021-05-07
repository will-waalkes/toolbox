# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:36:30 2017

@author: willwaalkes
"""

import os
import sys
import batman
import numpy as np
import matplotlib.pyplot as plt

#Next we create a TransitParams object to store the physical parameters describing the transit:

params = batman.TransitParams()
params.t0 = 0.                       #time of inferior conjunction
params.per = 1.628930                #orbital period
#params.rp = 0.0512                   #planet radius (in units of stellar radii)
params.rp = .1
params.a = 16.                       #semi-major axis (in units of stellar radii)
params.inc = 90.                     #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1956, 0.3700]          #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

#Note that for circular orbits, batman uses the convention params.w = 90. 
#The units for params.t0 and params.per can be anything as long as they are consistent.

#We also need to specify the times at which we wish to calculate the model:
c = 500
t = np.linspace(-0.05, 0.05, c)

#Using these parameters, we initialize the model and calculate a model light curve:

N = np.logspace(0,6,7)

m = batman.TransitModel(params, t)    #initializes model

#N = 1
#print np.random.normal(0,1,len(t))*np.sqrt(N)#*np.sqrt(N[i])/N[i]
#os.sys.exit()
for i in range(len(N)):
    flux = m.light_curve(params)*N[i]          #calculates light curve

    #Now simulate data points

#    sim = flux + np.random.normal(0,1,len(t))*np.sqrt(N[i])/N[i]
    sim = (np.random.poisson(flux,c))#/N[i]-1)

    plt.plot(t, flux)
    #plt.scatter(t, sim)
    plt.xlabel("Time from central transit")
    plt.ylabel("Relative flux")
    plt.title('Counts per minute (N):')
    plt.title(N[i], loc = 'right')
    plt.show()
    plt.clf()

os.sys.exit() 

plt.plot(N,np.sqrt(N))   
plt.xlabel("Number of Photons per exposure")
plt.xscale("log")
plt.ylabel("Expected Noise")
plt.yscale("log")
plt.title("Raw Noise")#plt.show()

plt.plot(N,np.sqrt(N)/N)   
plt.xlabel("Number of Photons per exposure")
plt.xscale("log")
plt.ylabel("Expected Noise")
plt.yscale("log")
plt.title("Fractional Noise")
plt.show()