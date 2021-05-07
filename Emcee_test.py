# -*- coding: utf-8 -*-
"""
Created on Thu May 18 15:46:26 2017

@author: willwaalkes
"""

import numpy as np
import emcee
import matplotlib.pyplot as pl

## We want to code something that is the probability - the inverse sum (sum^-1)
# cov = sum here. 

def lnprob(x, mu, icov):
    diff = x-mu
    return -np.dot(diff,np.dot(icov,diff))/2.0 #np.dot is the dot product
    
ndim = 50 #The number of dimensions....What is it for exoplanet data? 1D per parameter?

means = np.random.rand(ndim) #Constructing random means. Is my mean going to be the model?

cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov,cov)

icov = np.linalg.inv(cov)

nwalkers = 250
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[means, icov])

#'Burn in' the walkers so that they start scattered around
pos, prob, state = sampler.run_mcmc(p0, 100)
sampler.reset()

sampler.run_mcmc(pos, 1000)

for i in range(ndim):
    pl.figure()
    pl.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
    pl.title("Dimension {0:d}".format(i))

pl.show()

print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(sampler.acceptance_fraction)))