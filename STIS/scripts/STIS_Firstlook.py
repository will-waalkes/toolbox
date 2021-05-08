# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 12:05:06 2017

@author: willwaalkes
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import astropy.wcs as wcs
import stistools.x1d
import os
import sys

GI436_dir_inp = './Data/GI436b/2DSpect/'
GJ1132_dir_inp = './Data/GJ1132B/2DSpect/'

hdus = fits.open(GI436_dir_inp+'oc2403020_x2d.fits')

primary = hdus[0].data  # Primary (NULL) header data unit
img = hdus[1].data      # Intensity data
err = hdus[2].data      # Error per pixel
dq = hdus[3].data       # Data quality per pixel

plt.imshow(img, origin = 'lower',vmin = 0, vmax = 1.4e-12)
plt.colorbar()
plt.ylim(485,500)
plt.xlim(400,550)
print img.min()  # Call the min() method on the img object.
print img.max()


###### PLOT THE SPECTRUM AND COMPARE TO THE BACKGROUND

Spectrum = img[485:500,400:550].sum(axis=0)
plt.figure()
plt.plot(Spectrum)
plt.title('Spectrum')

background1 = img[510:525,400:550].sum(axis=0)
plt.figure()
plt.plot(background1)
plt.title('background1')

background2 = img[460:475,400:550].sum(axis=0)
plt.figure()
plt.plot(background2)
plt.title('background2')

background = (background1+background2)/2
plt.figure()
plt.plot(background)
plt.title('background')

calibrated = Spectrum-background
plt.figure()
plt.plot(calibrated)
plt.title('Background Subtracted')

###### NOW, LOAD A 1D SPECTRUM FROM THE STIS PIPELINE AND COMPARE

hdu = astropy.io.fits.open(GI436_dir_inp+'oc2403020_x1d.fits')
data = hdu[1].data
header = hdu[1].header
    
w = data['wavelength'].flatten()
f = data['flux'].flatten()
e = data['error'].flatten()

'''
os.path.basename('/Volumes/HUBBLEDATA/GI436b/1DSpect/oc2403040_flt.fits').replace('flt', 'x1d')
stistools.x1d.x1d('/Volumes/HUBBLEDATA/GI436b/1DSpect/oc2403040_flt.fits',
                  '/Volumes/HUBBLEDATA/GI436b/1DSpect/oc2403040_x1d.fits', maxsrch=0, extrsize=7,
                a2center=495, bk1offst=-20, bk2offst=20, backord=1,verbose=True)
'''

plt.figure()
#plt.plot(w,f)
plt.errorbar(w,f,yerr=e)
plt.xlabel('Wavelength (angstroms)')
plt.ylabel('Flux (erg/s/cm**2/angstrom)')
plt.xlim(1210,1220)
#plt.ylim(-0.5e-14,1e-14)
