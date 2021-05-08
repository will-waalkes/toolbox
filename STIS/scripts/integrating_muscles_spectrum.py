#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:45:04 2018

@author: willwaalkes
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

spectrum = fits.open('/Users/willwaalkes/Desktop/hlsp_muscles_multi_multi_gj1214_broadband_v22_var-res-sed.fits')

data = spectrum[1].data

wave_mid = data['WAVELENGTH']
wave_left = data['WAVELENGTH0']
wave_right = data['WAVELENGTH1']

flux = data['FLUX']

plt.plot(wave_mid,flux)
plt.xlim(0,1000)
plt.ylim(0,0.6e-14)

print wave_mid[4000]

dlam = wave_right-wave_left

r = 3.974e+19 #42 ly to cm
integral = np.sum(flux[0:4000]*dlam[0:4000])*4.0*np.pi*r**2.0
print integral*(14.6/12.04)**2.0
