#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:40:31 2017

@author: willwaalkes
"""

import os
from glob import glob
import astropy.io.fits
import matplotlib.pyplot as plt
import stistools.x1d 
import scipy as sc
from scipy import integrate
import numpy as np
from astropy.time import Time
from astropy import coordinates as coord, units as u

dir_in = '~/Research/GJ1132b/TIME-TAG/'

#hdu = astropy.io.fits.open(dir_in+'oda001010_tag.fits')
hdu = astropy.io.fits.open(dir_in+'Calibration_Files/oda001010_x2d.fits')
print hdu.info()
data = hdu[1].data
header = hdu[0].header
#GTI = hdu[2]
#print header
#print np.shape(data[0])
#print np.shape(hdu)
#print data
#print header
#plt.imshow(data)
#plt.colorbar()

#w = data['wavelength'].flatten()
#f = data['flux'].flatten()
#e = data['error'].flatten()
#StartTime = header['TEXPSTRT']
#EndTime = header['TEXPEND']
Duration = header['TEXPTIME']

#time = data['TIME']
#print time[1:]-time[:-1]

print Duration
#print 'DATA ',data
#print "SOMETHING ", GTI
#plt.scatter(data["AXIS1"],data["AXIS2"], s = 1, alpha=0.1)
#plt.show()
