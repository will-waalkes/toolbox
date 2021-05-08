#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:03:26 2017

@author: willwaalkes
"""

import os
from glob import glob
import astropy.io.fits
import matplotlib.pyplot as plt
import stistools.x1d 

'''
A re-extraction pipeline for STIS data of GI436b and GJ1132b.

This pipeline takes the _flt.fits 2-D spectral images and re-calibrates them into an x1d file
using stistools.

Why the _flt.fits file? These are uncorrected flat-fielded images, which will allow more
accurate data processing because they are unbinned.

At the end, this program should produce and plot a new set of 1D spectra for the input sources,
allowing the creation of transit light curves.
'''

dir_out = '~/Research/GJ1132b/Calibration_Output/'
dir_TT = '~/Research/GJ1132b/TIME-TAG/Calibration_Files/Separated_Exposures/'
dir_oref = '~/Research/GJ1132b/TIME-TAG/Calibration_Files/'
rest = 1215.67 #Rest wavelength in angstroms for these spectra

def reextract(source, filename, center, dir_inp, file_out):
    
    '''
    Uncomment the first part if a NEW calibration is necessary. The rest of the program plots the spectra.
    '''
    
    
    os.environ['oref'] = dir_oref
    print os.environ['oref']
    stistools.x1d.x1d(filename, file_out, maxsrch=4, extrsize=10,
            a2center=center, bk1offst=-20, bk2offst=20, backord=1,verbose=True)
    
    # Now plot the re-extracted spectrum
    hdu = astropy.io.fits.open(dir_out+source+'-1DSpect-'+str(i)+'.fits')
    data = hdu[1].data
    header = hdu[1].header
    
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()
    print 

    plt.figure()
    plt.plot(w,f)
    plt.errorbar(w,f,yerr=e)
    plt.xlabel('Wavelength (angstroms)')
    plt.ylabel('Flux (erg/s/cm**2/angstrom)')

    if source == 'GJ1132b':
        v = 35.0
        gj1132 =(1+v/3e5)*rest
        plt.axvline(gj1132, linestyle='--', alpha=0.3, color='gray', zorder=-100, label='GJ1132 rest-frame\n    ({0:.0f} km/s)'.format(v))
    else:
        plt.axvline(rest, linestyle='--',alpha=0.3, color='red', zorder=-100, label='heliocentric\n    ({0:.0f} km/s)'.format(0))
    plt.legend(frameon=False, loc='upper right', fontsize=10)
    plt.title(source)
    plt.title(i, loc = 'right')
    plt.xlim(1214,1218)
    plt.ylim(-.5e-14,)
    plt.show()
    
def TTextract(source, filename, center, dir_inp, file_out):
    
    '''
    Uncomment the first part if a NEW calibration is necessary. The rest of the program plots the spectra.
    '''
    
    
    os.environ['oref'] = dir_oref
    print os.environ['oref']
    stistools.x1d.x1d(filename, file_out, maxsrch=0, extrsize=10,
            a2center=center, bk1offst=-30, bk2offst=30, backord=1,verbose=True)
    
    #print file_out
    #os.sys.exit()
    
    # Now plot the re-extracted spectrum
    hdu = astropy.io.fits.open(dir_TT+source+'-TT1D-'+str(i)+'.fits')
    print dir_TT+source+'-TT1D-'+str(i)+'.fits'
    os.sys.exit()
    data = hdu[1].data
    header = hdu[1].header
    
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()

    plt.figure()
    plt.plot(w,f)
    plt.errorbar(w,f,yerr=e)
    plt.xlabel('Wavelength (angstroms)')
    plt.ylabel('Flux (erg/s/cm**2/angstrom)')

    v = 35.0
    gj1132 =(1+v/3e5)*rest
    plt.axvline(gj1132, linestyle='--', alpha=0.3, color='gray', zorder=-100, label='GJ1132 rest-frame\n    ({0:.0f} km/s)'.format(v))
    plt.legend(frameon=False, loc='upper right', fontsize=10)
    plt.title(source)
    plt.title(i, loc = 'right')
    plt.xlim(1214,1218)
    plt.ylim(-.5e-14,8e-14)
    plt.show()
    
'''
Call the calibration and plotting program with the specific path for your source
'''

GI436_x2d = glob('/Volumes/HUBBLEDATA/GI436b/'+'*0_flt.fits')
GJ1132_x2d = glob('/Volumes/HUBBLEDATA/GJ1132b/'+'*0_flt.fits')
GJ1132_TT = glob(dir_TT+'*_flt.fits')

i = 0

for filename in GI436_x2d: 
    i += 1
    #print i
    #reextract('GI436b', filename, 383.0, '/Volumes/HUBBLEDATA/GI436b/', dir_out+'GI436b-1DSpect-'+str(i)+'.fits')

i = 0
    
for filename in GJ1132_x2d:
    i += 1
    #print i
    #reextract('GJ1132b', filename, 405.0, '/Volumes/HUBBLEDATA/GI436b/', dir_out+'GJ1132b-1DSpect-'+str(i)+'.fits')

i = 0

for filename in GJ1132_TT:
    i += 1
    print i
    TTextract('GJ1132b', filename, 405.0, dir_oref, dir_TT+'GJ1132b-TT1D-'+str(i)+'.fits')
    #os.sys.exit()

#reextract('GJ1132b', dir_out+'oda001010_flt.fits', 405.0, '/Volumes/HUBBLEDATA/GJ1132b/', dir_out+'GJ1132b-1DSpect.fits')

