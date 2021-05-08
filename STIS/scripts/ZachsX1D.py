#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:46:02 2017

@author: willwaalkes
"""

import os, glob
import astropy.io.fits
import matplotlib.pyplot as plt
import stistools.x1d 

# identify files
#extractedpath = 'reextracted/'
#preface = '/Users/zkbt/Cosmos/Data/HST/GJ1132/STIS/'
#files = glob.glob(os.path.join(preface, 'oda0010*flt.fits'))

# set the extraction center
center = 493 #405 for GJ1132b? 493 for GI436b

# define an orbit
'''
def __init__(self, filename):
    self.fltfilename = filename
    self.reextract(filename)
    self.load()

def load(self): #loads an x1d file and extracts the wavelength, flux, error arrays
    hdu = astropy.io.fits.open(self.x1dfilename)
    self.data = hdu[1].data
    self.header = hdu[1].header
    
    self.w = self.data['wavelength'].flatten()
    self.f = self.data['flux'].flatten()
    self.e = self.data['error'].flatten()
    
def x1dfilename(self):
    return extractedpath + os.path.basename(self.fltfilename).replace('flt', 'x1d')
'''
    
def reextract(filename, a2center=405.0):
    os.environ['oref'] = '/Volumes/HUBBLEDATA/GJ1132b/'
    #
    print os.environ['oref']
    #reload(stistools.x1d)
    stistools.x1d.x1d(filename, 'test_x1d.fits', maxsrch=0, extrsize=7,
            a2center=a2center, bk1offst=-20, bk2offst=20, backord=1,
            verbose=True)#stistools.x1d
    #print os.environ['oref']
'''       
def plot(self):
    w, f, e = self.w, self.f, self.e
    d = self.data
    gs = plt.matplotlib.gridspec.GridSpec(2,1, hspace=0.2)
    axf = plt.subplot(gs[1])
    axf.errorbar(w,f, e, color='black', linewidth=0, elinewidth=2, capthick=2)
    axf.axhline(0, zorder=-100, color='gray', alpha=0.5)
    axf.get_xaxis().get_major_formatter().set_useOffset(False)
    axf.get_yaxis().get_major_formatter().set_useOffset(False)

    plt.xlabel('Wavelength (angstroms)')
    plt.ylabel('Flux (erg/s/cm**2/angstrom)')

    axc = plt.subplot(gs[0], sharex=axf)
    b = d['background'].flatten()
    g = d['gross'].flatten()
    n = d['net'].flatten()
    axc.plot(w,g, alpha=0.5, color='gray', label='total')
    axc.plot(w, b, alpha=0.3, color='red', label='background')
    axc.plot(w, n, alpha=1.0, color='black', label='subtracted')
    axc.legend(frameon=False, loc='upper left', fontsize=10)
    plt.ylabel('Flux (counts/s)')
    plt.setp(axc.get_xticklabels(), visible=False)
    plt.xlim([1214,1218])

    rest = 1215.7
    v = 35.0
    gj1132 =(1+v/3e5)*rest
    axf.axvline(gj1132, linestyle='--', alpha=0.3, color='gray', zorder=-100, label='GJ1132 rest-frame\n    ({0:.0f} km/s)'.format(v))
    axf.axvline(rest, linestyle='--',alpha=0.3, color='red', zorder=-100, label='heliocentric\n    ({0:.0f} km/s)'.format(0))
    axc.axvline(gj1132, linestyle='--',alpha=0.3, color='gray', zorder=-100)
    axc.axvline(rest, linestyle='--',alpha=0.3, color='red', zorder=-100)
    axf.legend(frameon=False, loc='upper right', fontsize=10)
    axf.set_ylim(-0.5e-14, 4e-14)

    axc.set_ylim(-0.005, 0.04)
    axc.set_title(self.x1dfilename)
        
    plt.savefig(self.x1dfilename.replace('fits', 'pdf'))
    plt.draw()
        
for f in files:
    s = stisOrbit(f)
    s.plot()
    plt.show()
    
    '''
    
reextract('oda004010_flt.fits')