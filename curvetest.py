#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 11:08:47 2017

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

'''
This is a test of overlaying a transit light curve onto my data.
'''

dir_in = '/Users/willwaalkes/Desktop/Colorado/ColoradoResearch/GJ1132b/Calibration_Output/'

def FluxInt(source, filename, rest, ra, dec, BJD_t0, period):
    
    '''
    FluxInt will convert UTC header times to BJD, convert BJD to orbital phase, and then integrate
    the red-shifted and blue-shifted wings of the Lyman-alpha spectra. Times and integrated values
    are stored to arrays for use in the next program.
    
    First, we must calculate the midpoint (between beginning and
    end of exposure) time, so that we can confine a long exposure to a single point.
    '''

    hdu = astropy.io.fits.open(dir_in+source+'-1DSpect-'+str(i)+'.fits')
    data = hdu[1].data
    header = hdu[1].header
    
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()
    UTC_date = header['DATE-OBS'] #UTC date of exposure
    start_time = header['TIME-OBS'] #UTC Time of exposure start
    EXP_time = header['EXPTIME'] #Exposure time in seconds
    
    exposuretimes = (EXP_time/(60.0*60.0))/2
    exposure.append(exposuretimes)
    
    f = f*1e14
    error = e*1e14
    error_sq = error**2 #Used for error propagation
        
    apo_location = [-105.820417*u.deg,32.780361*u.deg] #Using the APO coords since HST has none?
    exp_time = EXP_time/(3600.*24.) #put seconds into day
    t_start = Time(UTC_date+'T'+start_time, format = 'isot', scale='utc')
    start_time = t_start.jd #convert to jd times
    end_time = start_time + exp_time #in jd
    avg_time = (start_time+end_time)/2. #in jd
    avg_time_t = Time(avg_time, format='jd', scale='utc',location=(apo_location[0],apo_location[1]))  
     
    star_loc = coord.SkyCoord(ra,dec,unit=(u.hourangle,u.deg),frame='icrs')
    light_travel = avg_time_t.light_travel_time(star_loc)
    
    BJD_time = light_travel+avg_time_t.tdb
    BJD_time = BJD_time.jd
    
    '''
    Calculate the phase of each observation, based on linear relationship between BJD and period.
    '''
    
    Norbits = (BJD_time - BJD_t0)/period
    
    orbit = str(Norbits-int(Norbits))[1:]
    orbit = float(orbit)
    if orbit >= 0.5:
        phase = -1.0 + orbit
    else:
        phase = orbit
        
    transit_time = phase*(period*24.0) #This calculates how many hours before and after transit
    #print transit_time
    
    '''
    Now integrate the two lyman alpha wings
    '''   
    
    mid = np.amin(abs(w-rest))
    for j in range(len(w)):
        if abs(w[j] - rest) == mid:
            mid_pix=j
    lower_pix = mid_pix-12
    upper_pix = mid_pix+17
    
    blue = sc.integrate.trapz(f[lower_pix:(mid_pix+1)],dx=1)
    red = sc.integrate.trapz(f[mid_pix:(upper_pix+1)],dx=1)
    blue_err = np.sqrt(sc.integrate.trapz(error_sq[lower_pix:(mid_pix+1)],dx=1))
    red_err = np.sqrt(sc.integrate.trapz(error_sq[mid_pix:(upper_pix+1)],dx=1))
    
    blue_shift.append(blue)
    red_shift.append(red)
    blue_error.append(blue_err)
    red_error.append(red_err)
    phase_vals.append(phase)
    transit_times.append(transit_time)
        
def LightCurve(source, phase_vals, blue_shift, red_shift, blue_error, red_error):
    
    '''
    Plot the light curves
    '''
    
    plt.figure()
    if source == 'GI436b':
        #plt.scatter(phase_vals[0:4],blue_shift[0:4], label = 'visit 2')
        #plt.scatter(phase_vals[4:8],blue_shift[4:8], label = 'visit 3')
        #plt.errorbar(phase_vals[0:4],blue_shift[0:4], yerr = blue_error[0:4], fmt='o')
        #plt.errorbar(phase_vals[4:8],blue_shift[4:8], yerr = blue_error[4:8], fmt='o')
        plt.axvline(-0.00788, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.axvline(0.00788, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.ylim(15,80)
    if source == 'GJ1132b':
        #plt.scatter(phase_vals[0:5],blue_shift[0:5], label = 'visit 1')
        #plt.scatter(phase_vals[7:12],blue_shift[7:12], label = 'visit 2')
        #plt.errorbar(phase_vals[0:5],blue_shift[0:5], yerr = blue_error[0:5], fmt='o')
        #plt.errorbar(phase_vals[7:12],blue_shift[7:12], yerr = blue_error[7:12], fmt='o')
        plt.axvline(-0.02, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.axvline(0.02, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.ylim(-2.0,28.0)
    
    plt.scatter(phase_vals,blue_shift, color = 'k')
    plt.errorbar(phase_vals, blue_shift, yerr = blue_error, color = 'k', fmt = 'o')
    plt.legend(frameon=False, loc='upper left', fontsize=10)
    plt.title('Ly-alpha Blue Wing')
    plt.title(source, loc='left')
    plt.xlabel('Phase')
    plt.xlabel('Time since mid-transit (hr)')
    plt.ylabel('Flux ($x10^{-14}$ erg/s/cm$^2$/angstrom)')
    plt.show()
    
    plt.figure()
    if source == 'GI436b':
        #plt.scatter(phase_vals[0:4],red_shift[0:4], label = 'visit 2')
        #plt.scatter(phase_vals[4:8],red_shift[4:8], label = 'visit 3')
        #plt.errorbar(phase_vals[0:4],red_shift[0:4], yerr = red_error[0:4], fmt='o')
        #plt.errorbar(phase_vals[4:8],red_shift[4:8], yerr = red_error[4:8], fmt='o')
        plt.axvline(-0.00788, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.axvline(0.00788, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.ylim(15,80)
    if source == 'GJ1132b':
        #plt.scatter(phase_vals[0:5],red_shift[0:5], label = 'visit 1')
        #plt.scatter(phase_vals[7:12],red_shift[7:12], label = 'visit 2')
        #plt.errorbar(phase_vals[0:5],red_shift[0:5], yerr = red_error[0:5], fmt='o')
        #plt.errorbar(phase_vals[7:12],red_shift[7:12], yerr = red_error[7:12], fmt='o')
        plt.axvline(-0.02, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.axvline(0.02, linestyle='--', alpha=1.0, color='black', zorder=-100)
        plt.ylim(-2.0,28.0)
        
    plt.scatter(phase_vals,red_shift, color = 'k')
    plt.errorbar(phase_vals,red_shift, yerr = red_error, color = 'k', fmt='o')
    plt.legend(frameon=False, loc='lower left', fontsize=10)
    plt.title('Ly-alpha Red Wing')
    plt.title(source, loc='left')
    #plt.xlabel('Phase')
    plt.xlabel('Time since mid-transit (hr)')
    plt.ylabel('Flux ($x10^{-14}$ erg/s/cm$^2$/angstrom)')
    plt.show()

'''
CALL WHICH PROCESS YOU WANT TO RUN
'''

GI436_x1d = glob(dir_in+'*GI436*.fits')
GJ1132_x1d = glob(dir_in+'*GJ1132*.fits')

blue_shift = []
red_shift = []
blue_error = []
red_error = []
phase_vals = []
transit_times = []
exposure = []
i = 0

for filename in GI436_x1d: 
    i += 1
    #SpectImage('GI436b', filename, 383.0, 1215.67)
    FluxInt('GI436b', filename, 1215.67, '11:42:11.0', '+26:42:23', 2454510.80162, 2.643850)

LightCurve('GI436b', phase_vals, blue_shift, red_shift, blue_error, red_error)

blue_shift = []
red_shift = []
blue_error = []
red_error = []
phase_vals = []
transit_times = []
exposure = []
i = 0
    
for filename in GJ1132_x1d:
    i += 1
    #SpectImage('GJ1132b', filename, 405.0, 1215.8118)
    FluxInt('GJ1132b', filename, 1215.8118, '10:14:51.77', '-47:09:24.1', 2457184.55786, 1.62893)

LightCurve('GJ1132b', phase_vals, blue_shift, red_shift, blue_error, red_error)