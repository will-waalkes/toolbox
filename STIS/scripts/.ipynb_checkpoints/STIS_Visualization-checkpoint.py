#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 16:18:56 2017

@author: willwaalkes
"""
import os
import sys
from glob import glob
import astropy.io.fits
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits.hdu.hdulist as h
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii

'''
This script will visualize 1d spectra created from the STIS pipeline from x2d files. 
In particular, this will allow the visualization of the lyman-alpha spectrum of a star,
showing both individual time-series spectra and a composite stacked spectrum, with error propagation.
'''

c = 3.0e5 #km/s
dir_in = '~/Research/GJ1132b/Calibration_Output/'
desktop = '~/Desktop/'

def SingleSpectImage(source, filename, center, rest, v_rest):
    
    '''
    SpectImage reads in the x1d files and plots the Ly-alpha spectrum. This includes a calculation
    for where (pixel values) the blue-shift and red-shifted integrals will be calculated.
    '''

    hdu = astropy.io.fits.open(dir_in+source+'-1DSpect-'+str(i)+'.fits')
    data = hdu[1].data
    header = hdu[1].header
    
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()
    UTCDate = header['DATE-OBS']
    UTCTime = header['TIME-OBS']
    
    f = f*1e14
    e = e*1e14
    
    Composite_f = np.loadtxt('Composite_f.txt')
    Composite_err = np.loadtxt('Composite_err.txt')
    
    v_space = c*(w-rest)/w
    
    V_geocor_low = c*(1215.58-rest)/1215.58
    V_geocor_high = c*(1215.81-rest)/1215.75
    
    '''
    Here I am choosing integration regions based on Ehrenreich+ 2015. Blue region is -120 -> -40km/s.
    The red region is 30 -> 200 km/s. Should it be the same for GJ1132b?
    
    update - have now changed integraton zones to go down to 0 km/s
    '''
       
    blue_upper = np.amin(abs(v_space))
    for j in range(len(v_space)):
        if abs(v_space[j]) == blue_upper:
            upper_blue_pix = j
    blue_lower = np.amin(abs(v_space+200.0))
    for j in range(len(v_space)):
        if abs(v_space[j] + 200.0) == blue_lower:
            lower_blue_pix = j
    red_upper = np.amin(abs(v_space - 200.0))
    for j in range(len(v_space)):
        if abs(v_space[j] - 200.0) == red_upper:
            upper_red_pix = j
    red_lower = np.amin(abs(v_space))
    for j in range(len(v_space)):
        if abs(v_space[j]) == red_lower:
            lower_red_pix = j
       
    plt.figure()
    plt.plot(v_space,f)
    plt.errorbar(v_space,f,yerr=e, color = 'k')
    plt.xlabel('Velocity (km s$^{-1}$)')
    plt.ylabel('Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$ $\AA$)')
    
    plt.errorbar(v_space,Composite_f,yerr=Composite_err,
                 zorder=-1000,alpha=0.8,color='red',linestyle='--',label='Average Spectrum')

    plt.axvspan(V_geocor_low,V_geocor_high,label='Geocoronal Emission',color='orange',alpha=0.1,zorder=-1000) 
    if source == 'GJ1132b':
        plt.ylim(-1.0,5)
    if source == 'GI436b':
        plt.ylim(-2.0,12.0)
        
    plt.legend(frameon=False, loc='upper right', fontsize=10)
    plt.title(source)
    plt.title(UTCDate+' '+UTCTime, loc = 'right')
    plt.xlim(-220.0,220.0)
    plt.savefig(fname=desktop+source+'-'+str(UTCDate)+str(UTCTime)+'-spect.png',clobber=True) #Saving images to turn into a gif
    plt.show()
    
def CompositeImage(source, filename, center, rest, v_rest):
    
    
    '''
    CompositeImage will create a single plot with each exposure's spectrum, and then a
    final averaged spectrum.
    '''

    hdu = astropy.io.fits.open(dir_in+source+'-1DSpect-'+str(i)+'.fits')
    data = hdu[1].data
    header = hdu[1].header
    
    print(dir_in+source+'-1DSpect-'+str(i)+'.fits')
    print('-----------')
    
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()
    e1 = data['net_error'].flatten()
    phots = data['net'].flatten()
    phot_for_err = data['gross'].flatten()
    e2 = np.sqrt(np.abs(phot_for_err)*2095.0)/2095.0
    back = data['background'].flatten()
    
    '''
    plt.figure()
    plt.plot(w,back)
    plt.plot(w,phot_for_err)
    plt.axhline(0,color='purple')
    plt.title('background vs wavelength')
    plt.xlim(1215,1217)
    plt.ylim(-0.002,0.005)
    plt.show()
    '''
    
    #plt.figure()
    #plt.scatter(e2*2095.0,phot_for_err*2095.0)
    #plt.xlabel('Net Error')
    #plt.ylabel('sqrt(N/s) Error')
    #plt.title('error vs gross phot')
    #plt.show()
    
    #plt.figure()
    #plt.hist((e2-e1),bins=30)
    
    #plt.figure()
    #plt.hist(e2,bins=30)
    #print e
    #phots = data['gross'].flatten()
    
    #e = e2
    
    '''A note about the errors: e2 is the proper way to do the photon counting errors. e1 is the pipeline
    calculation of this error. e is the flux error from the pipeline'''
        
    f = f*1e14
    #f = phots
    f_array.append(f)
    e = e*1e14
    error_array.append(e**2)
    inv_err.append(1/(e**2))
    
    v_space = c*(w-rest)/w
    
    V_geocor_low = c*(1215.58-rest)/1215.58
    V_geocor_high = c*(1215.825-rest)/1215.75
    
            
    '''
    Here I am choosing integration regions based on Ehrenreich+ 2015. Blue region is -120 -> -40km/s.
    The red region is 30 -> 200 km/s. Should it be the same for GJ1132b?
    '''
       
    blue_upper = np.amin(abs(v_space))
    for j in range(len(v_space)):
        if abs(v_space[j]) == blue_upper:
            upper_blue_pix = j
    blue_lower = np.amin(abs(v_space+160.0))
    for j in range(len(v_space)):
        if abs(v_space[j] + 160.0) == blue_lower:
            lower_blue_pix = j
    red_upper = np.amin(abs(v_space - 160.0))
    for j in range(len(v_space)):
        if abs(v_space[j] - 160.0) == red_upper:
            upper_red_pix = j
    red_lower = np.amin(abs(v_space))
    for j in range(len(v_space)):
        if abs(v_space[j]) == red_lower:
            lower_red_pix = j 
            
    '''
    This plotting function will overplot each individual spectrum
    '''
    
        
    '''    
    plt.figure()
    for j in range(len(f_array)):
        if j == 3:
            plt.plot(v_space,f_array[j],color='red',label='In transit',zorder=100)
            plt.errorbar(v_space,f_array[j],yerr=error_array[j],color='red',zorder=100,fmt='--')
        #if j == 10:
        #    plt.plot(v_space,f_array[j],color='cyan',label='Visit 2 in transit',zorder=100)
        #    plt.errorbar(v_space,f_array[j],yerr=error_array[j],color='cyan',zorder=100,fmt='--')
            
        if j == 0:
            plt.plot(v_space,f_array[j],color='gray',label='Out of transit')
            plt.errorbar(v_space,f_array[j],yerr=error_array[j],color='gray',fmt='--')
            plt.axvspan(V_geocor_low,V_geocor_high,label='Geocoronal Emission',color='blue',alpha=0.15,zorder=-1000)        
        
        if (j == 1 or j == 2 or j ==4 or j == 5 or j == 6):
            plt.plot(v_space,f_array[j],color = 'gray')
            plt.errorbar(v_space,f_array[j],yerr=error_array[j],color='gray',fmt='--')
            
        #if j == 7:
        #    plt.plot(v_space,f_array[j],color='black',label='Visit 2 out of transit')
        #    plt.errorbar(v_space,f_array[j],yerr=error_array[j],color='black',fmt='--')
        #if (j == 8 or j == 9 or j == 11 or j == 12 or j == 13):
        #    plt.plot(v_space,f_array[j],color = 'black')
        #    plt.errorbar(v_space,f_array[j],yerr=error_array[j],color='black',fmt='--')
            
    plt.xlabel('Velocity (relative to GJ1132) (km s$^{-1}$)',fontsize=15)
    plt.ylabel('Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$ $\AA$)',fontsize=15)
    #plt.ylabel('photons/s')

    #plt.axvline(0, linestyle='--', alpha=1.0, color='orange', zorder=-100,label='Rest-frame')
    #plt.axvline(v_space[lower_blue_pix], linestyle='--', alpha=0.3, color='blue', zorder=-100)
    #plt.axvline(v_space[upper_red_pix], linestyle='--', alpha=0.3, color='red', zorder=-100)
    #plt.axvline(v_space[upper_blue_pix], linestyle='--', alpha=0.3, color='blue', zorder=-100)
    #plt.axvline(v_space[lower_red_pix], linestyle='--', alpha=0.3, color='red', zorder=-100)
    plt.axhline(0.0,alpha=0.5,color='purple',zorder=100)
    
    #if source == 'GJ1132b':
    #    plt.ylim(-1.0,5.5)
        
    plt.legend(frameon=False, loc='upper left', fontsize=10)
    plt.title('Visit 1 Spectra',fontsize=18)
    plt.xlim(-220.0,220.0)
    plt.tight_layout()
    plt.savefig(desktop+'V1_multispectra.pdf',clobber=True)
    plt.show()
    '''
    
    
    '''
    Here I will create the total average spectrum and weighted errors.
    '''
    
    
    #Lya_model = np.loadtxt('/Users/willwaalkes/Desktop/Lya_model_best.txt')
    #Lya_model = Lya_model*1e14
    #wave_to_fit = np.loadtxt('/Users/willwaalkes/Desktop/model_waves.txt')
    #wave_to_vspace = c*(wave_to_fit-rest)/wave_to_fit
    
    Composite_f = np.sum(np.array(f_array)*np.array(inv_err),axis=0)/np.sum(inv_err,axis=0)
    Composite_err = np.sqrt((1/np.sum(inv_err,axis=0)))
    
    
    data = Table([Composite_f, w, Composite_err,], names=['flux', 'wave', 'error'])
    ascii.write(data, 'GJ1132b_avg_spect.dat',overwrite=True)
    
    plt.figure()
    plt.errorbar(v_space,Composite_f,yerr=Composite_err, color = 'k')
    #plt.plot(wave_to_vspace,Lya_model,zorder=-1000,alpha=0.5,color='gray',linestyle='--')
    
    plt.xlabel('Velocity (relative to GJ1132) (km s$^{-1}$)',fontsize=15)
    plt.ylabel('Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$ $\AA$)',fontsize=15)
    #plt.ylabel('Photons/s')

    #plt.axvline(0, linestyle='--', alpha=1.0, color='gray', zorder=-100,label='GJ1132 Rest Frame')
    plt.axvline(v_space[lower_blue_pix], linestyle='--', alpha=0.3, color='blue', zorder=-100)
    plt.axvline(v_space[upper_red_pix], linestyle='--', alpha=0.3, color='red', zorder=-100)
    #plt.axvline(v_space[upper_blue_pix], linestyle='--', alpha=0.3, color='blue', zorder=-100)
    #plt.axvline(v_space[lower_red_pix], linestyle='--', alpha=0.3, color='red', zorder=-100)
    #plt.ylim(-1.0,)
        
    plt.legend(frameon=False, loc='upper right', fontsize=10)
    plt.title('Stacked Spectrum',fontsize=18)
    plt.xlim(-220.0,220.0)
    
    plt.fill_between(v_space[lower_blue_pix:upper_blue_pix-3],
                     Composite_f[lower_blue_pix:upper_blue_pix-3],
                     color='royalblue',alpha=0.8)
    plt.fill_between(v_space[lower_red_pix+1:upper_red_pix+1],
                     Composite_f[lower_red_pix+1:upper_red_pix+1],
                     color='crimson',alpha=0.6)
    
    #print('Blue wave range:',v_space[lower_blue_pix:upper_blue_pix-3])
    #print('Red wave range:',v_space[lower_red_pix+1:upper_red_pix+1])
    
    plt.tight_layout()
    plt.savefig(desktop+'Cumulative_spectrum.pdf', clobber=True)
    plt.show()
    
    #if i == 7:
    #    os.sys.exit()
    
    
    
#    mask1 = w >= 1214
#    w = w[mask1]
#    Comp_f = Composite_f[mask1]
#    Comp_err = Composite_err[mask1]
#    
#    mask2 = w <= 1218
#    w = w[mask2]
#    Comp_f = Comp_f[mask2]/1.0e14
#    Comp_err = Comp_err[mask2]/1.0e14
#    
#    
#    data = Table([w,Comp_f,Comp_err],names=['wave','flux (photons/s)','error (photons/s)'])
#    plt.figure()
#    plt.errorbar(data['wave'], data['flux'],yerr=data['error'])
#    plt.title('Final Spectrum')
#    plt.show()
#    ascii.write(data, '/Users/willwaalkes/Desktop/GJ1132b_photdata.txt',overwrite=True)
#    print Comp_err
    
    
'''
WRITE THE COMPOSITE FLUX AND ERROR TO A NEW HDULIST TO SAVE AS AN X1D
'''    
    
    
'''
CALL WHICH PROCESS YOU WANT TO RUN
'''

GJ1132_x1d = glob(dir_in+'*GJ1132*.fits')

blue_shift = []
red_shift = []
f_array = []
error_array = []
inv_err = []
i = 0
    
for filename in GJ1132_x1d:
    i += 1
    #SingleSpectImage('GJ1132b', filename, 405.0, 1215.8118, 35.0)
    CompositeImage('GJ1132b', filename, 405.0, 1215.8118, 35.0)