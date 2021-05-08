"""
Created on Thu Jun 22 09:38:53 2017

@author: willwaalkes
"""

import os
from glob import glob
import astropy.io.fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy import coordinates as coord, units as u

'''
Spect_processing is used to visualize HST STIS x1d spectral data at the Lyman-alpha wavelength.
This script contains programs that will plot the spectrum, integrate the blue-shifted and
red-shifted wings, and generate light curves. Also included are calculations for converting UTC
to BJD and determining orbital phase.

For GI 436b, the rest wavelength is taken to be heliocentric - 1215.67 Angstrom
For GJ 1132b, the rest wavelength is calculated from GJ_rest = (1+v/3e5)*rest where v is 35 km/s
'''

c = 3.0e5 #km/s
dir_in = '~/Research/GJ1132b/Calibration_Output/'
dir_TT = '~/Research/GJ1132b/TIME-TAG/Calibration_Files/Separated_Exposures/'

def FluxInt(source, filename, rest, ra, dec, BJD_t0, period, v_rest):
    
    '''
    FluxInt will convert UTC header times to BJD, convert BJD to orbital phase, and then integrate
    the red-shifted and blue-shifted wings of the Lyman-alpha spectra. Times and integrated values
    are stored to arrays for use in the next program.
    
    First, we must calculate the midpoint (between beginning and
    end of exposure) time, so that we can confine a long exposure to a single point.
    '''

    hdu = astropy.io.fits.open(filename)
    data = hdu[1].data
    header = hdu[1].header
    
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()
    UTC_date = header['DATE-OBS'] #UTC date of exposure
    start_time = header['TIME-OBS'] #UTC Time of exposure start
    EXP_time = header['EXPTIME'] #Exposure time in seconds
    Fullexp_start = header['EXPSTART']
    
    f = f*1.0e14
    error = e*1.0e14
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
    Calculate the phase of each observation, based on the linear relationship between BJD and period.
    '''
    
    Norbits = (BJD_time - BJD_t0)/period
    
    orbit = str(Norbits-int(Norbits))[1:]
    orbit = float(orbit)
    if orbit >= 0.5:
        phase = -1.0 + orbit
    else:
        phase = orbit
        
    transit_time = phase*(period*24.0) #This calculates how many hours before and after transit
    
    '''
    Now integrate the two lyman alpha wings
    '''   
    
    v_space = c*(w-rest)/w
            
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
    Integrate the flux over the blue-shifted and red-shifted wavelength ranges.
    '''
    
    blue = np.sum(f[lower_blue_pix:upper_blue_pix-3]*0.0533)
    blue_err = np.sqrt(np.sum(error_sq[lower_blue_pix:upper_blue_pix-3]))*0.0533
    red = np.sum(f[lower_red_pix+1:upper_red_pix+2]*0.0533)
    red_err = np.sqrt(np.sum(error_sq[lower_red_pix+1:upper_red_pix+2]))*.0533
        
    if float(Fullexp_start) <= 57800.0:
        transit_times_v1.append(transit_time)
        blue_error_v1.append(blue_err)
        red_error_v1.append(red_err)
        red_shift_v1.append(red)
        blue_shift_v1.append(blue)

        dat_b = Table([blue_shift_v1,blue_error_v1,transit_times_v1],
                      names=['Flux','Error','Time_hr'])
        dat_r = Table([red_shift_v1,red_error_v1,transit_times_v1],
                      names=['Flux','Error','Time_hr'])

        ascii.write(dat_b,'full_visit1_blue.txt',overwrite=True)
        ascii.write(dat_r,'full_visit1_red.txt',overwrite=True)

    else:
        transit_times_v2.append(transit_time)
        blue_error_v2.append(blue_err)
        red_error_v2.append(red_err)
        red_shift_v2.append(red)
        blue_shift_v2.append(blue)

        dat_b = Table([blue_shift_v2,blue_error_v2,transit_times_v2],
                      names=['Flux','Error','Time_hr'])
        dat_r = Table([red_shift_v2,red_error_v2,transit_times_v2],
                      names=['Flux','Error','Time_hr'])

        ascii.write(dat_b,'full_visit2_blue.txt',overwrite=True)
        ascii.write(dat_r,'full_visit2_red.txt',overwrite=True)
    
def TT_FluxInt(source, filename, rest, ra, dec, BJD_t0, period, v_rest):
    
    '''
    FluxInt will convert UTC header times to BJD, convert BJD to orbital phase,
    and then integrate the red-shifted and blue-shifted wings of the Lyman-alpha spectra.
    Times and integrated values are stored to arrays for use in the next program.
    
    First, we must calculate the midpoint (between beginning and
    end of exposure) time, so that we can confine a long exposure to a single point.
    '''

    hdu = astropy.io.fits.open(filename)
    data = hdu[1].data
    header = hdu[1].header
    w = data['wavelength'].flatten()
    f = data['flux'].flatten()
    e = data['error'].flatten()
    TT_start = header['EXPSTART'] #In MJD
    TT_end = header['EXPEND'] #MJD
        
    f = f*1.0e14
    error = e*1.0e14
    error_sq = error**2 #Used for error propagation
        
    apo_location = [-105.820417*u.deg,32.780361*u.deg] #Using the APO coords since HST has none?
    start_time = TT_start + 2400000.5 #converting to JD
    end_time = TT_end + 2400000.5
    avg_time = (start_time+end_time)/2.
    avg_time_t = Time(avg_time, format='jd', scale='utc',location=(apo_location[0],apo_location[1]))  
     
    star_loc = coord.SkyCoord(ra,dec,unit=(u.hourangle,u.deg),frame='icrs')
    light_travel = avg_time_t.light_travel_time(star_loc)
    
    BJD_time = light_travel+avg_time_t.tdb
    BJD_time = BJD_time.jd
    
    '''
    Calculate the phase of each observation, based on the linear relationship
    between BJD and period.
    '''
    
    Norbits = (BJD_time - BJD_t0)/period
    
    orbit = str(Norbits-int(Norbits))[1:]
    orbit = float(orbit)
    if orbit >= 0.5:
        phase = -1.0 + orbit
    else:
        phase = orbit
        

    transit_time = phase*(period*24.0) #This calculates how many hours before and after transit
    
    '''
    Now integrate the two lyman alpha wings
    '''   
    
    v_space = c*(w-rest)/w
            
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
    Integrate the flux over the blue-shifted and red-shifted wavelength ranges.
    '''
    
    blue = np.sum(f[lower_blue_pix:upper_blue_pix-3]*0.0533)
    blue_err = np.sqrt(np.sum(error_sq[lower_blue_pix:upper_blue_pix-3]))*0.0533
    red = np.sum(f[lower_red_pix+1:upper_red_pix+2]*0.0533)
    red_err = np.sqrt(np.sum(error_sq[lower_red_pix+1:upper_red_pix+2]))*.0533
    
    if float(TT_start) <= 57800.0:
        transit_times_v1.append(transit_time)
        blue_error_v1.append(blue_err)
        red_error_v1.append(red_err)
        red_shift_v1.append(red)
        blue_shift_v1.append(blue)

        dat_b = Table([blue_shift_v1,blue_error_v1,transit_times_v1],
                      names=['Flux','Error','Time_hr'])
        dat_r = Table([red_shift_v1,red_error_v1,transit_times_v1],
                      names=['Flux','Error','Time_hr'])

        ascii.write(dat_b,'erg_visit1_blue.txt',overwrite=True)
        ascii.write(dat_r,'erg_visit1_red.txt',overwrite=True)

    else:
        transit_times_v2.append(transit_time)
        blue_error_v2.append(blue_err)
        red_error_v2.append(red_err)
        red_shift_v2.append(red)
        blue_shift_v2.append(blue)

        dat_b = Table([blue_shift_v2,blue_error_v2,transit_times_v2],
                      names=['Flux','Error','Time_hr'])
        dat_r = Table([red_shift_v2,red_error_v2,transit_times_v2],
                      names=['Flux','Error','Time_hr'])

        ascii.write(dat_b,'erg_visit2_blue.txt',overwrite=True)
        ascii.write(dat_r,'erg_visit2_red.txt',overwrite=True)
    

'''
CALL WHICH PROCESS YOU WANT TO RUN
'''

GI436_x1d = glob(dir_in+'*GI436*.fits')
GJ1132_x1d = glob(dir_in+'*GJ1132*.fits')
GJ1132_TT = glob(dir_TT+'*TT1D*.fits')

blue_shift = []
red_shift = []
blue_error = []
red_error = []
phase_vals = []
transit_times = []
i = 0

for filename in GI436_x1d: 
    i += 1
    #FluxInt('GI436b', filename, 1215.67, '11:42:11.0', '+26:42:23', 2454865.083208, 2.6438979, 9.6)

blue_shift_v1 = []
blue_shift_v2 = []
red_shift_v1 = []
red_shift_v2 = []
blue_error_v1 = []
blue_error_v2 = []
red_error_v1 = []
red_error_v2 = []
phase_vals_v1 = []
phase_vals_v2 = []
transit_times_v1 = []
transit_times_v2 = []
i = 0
    
for filename in GJ1132_x1d:
    i += 1
    #FluxInt('GJ1132b', filename, 1215.8118, '10:14:51.77', '-47:09:24.1', 2457184.55786, 1.62893, 35.0)

blue_shift_v1 = []
blue_shift_v2 = []
red_shift_v1 = []
red_shift_v2 = []
blue_error_v1 = []
blue_error_v2 = []
red_error_v1 = []
red_error_v2 = []
phase_vals_v1 = []
phase_vals_v2 = []
transit_times_v1 = []
transit_times_v2 = []
i = 0
    
for filename in GJ1132_TT:
    i += 1
    TT_FluxInt('GJ1132b', filename, 1215.8118, '10:14:51.77', '-47:09:24.1', 2457184.55786, 1.62893, 35.0)