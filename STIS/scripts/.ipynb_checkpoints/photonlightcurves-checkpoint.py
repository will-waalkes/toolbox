import os
from glob import glob
import astropy.io.fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy import coordinates as coord, units as u
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import batman
import scipy.stats as stats

c = 3.0e5 #km/s
dir_in = '~/Research/GJ1132b/Calibration_Output/'
dir_TT = '~/Research/GJ1132b/TIME-TAG/Calibration_Files/Separated_Exposures/'
desktop = '~/Desktop/'

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

    TT_start = header['EXPSTART'] #In MJD
    TT_end = header['EXPEND'] #MJD

    #exptime = (TT_end-TT_start)*24.0*60.0*60.0
    exptime = header['EXPTIME']

    flux = data['flux'].flatten()
    w = data['wavelength'].flatten()
    phot_rate = data['net'].flatten() #net counts per second
    phots_net = phot_rate*exptime #total net photons
    gross = data['gross'].flatten() #gross counts per second

    phots_gross = gross*exptime #total gross photons
    e = np.sqrt(np.abs(phots_gross))
    
    plt.scatter(w,flux/phot_rate,s=4,alpha=0.5)
    plt.xlabel('wavelength')
    plt.ylabel('flux/phot_rate')
    plt.xlim(1214,1217)
    plt.ylim(-2e-12,0.4e-11)

    error_sq = e**2 #Used for error propagation

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

    V_geocor_low = c*(1215.58-rest)/1215.58
    V_geocor_high = c*(1215.75-rest)/1215.75

    v_space = c*(w-rest)/w

    lower_blue_pix = 384
    upper_blue_pix = 392
    lower_red_pix = 397
    upper_red_pix = 410

    '''
    Integrate the flux over the blue-shifted and red-shifted wavelength ranges.
    '''

    blue_err = np.sqrt(np.sum(error_sq[lower_blue_pix:(upper_blue_pix+1)]))
    blue_gross = np.sum(phots_gross[lower_blue_pix:(upper_blue_pix+1)])
    blue_net = np.sum(phots_net[lower_blue_pix:(upper_blue_pix+1)])

    red_err = np.sqrt(np.sum(error_sq[lower_red_pix:(upper_red_pix+1)]))
    red_gross = np.sum(phots_gross[lower_red_pix:(upper_red_pix+1)])
    red_net = np.sum(phots_net[lower_red_pix:(upper_red_pix+1)])

    if float(TT_start) <= 57800.0:
        exptime_v1.append(exptime)
        transit_times_v1.append(transit_time)
        blue_error_v1.append(blue_err)
        red_error_v1.append(red_err)
        red_net_v1.append(red_net)
        blue_net_v1.append(blue_net)
        blue_gross_v1.append(blue_gross)
        red_gross_v1.append(red_gross)
        red_sky_v1 = np.array(red_gross_v1)-np.array(red_net_v1)
        blue_sky_v1 = np.array(blue_gross_v1)-np.array(blue_net_v1)

        dat_b = Table([blue_net_v1,blue_gross_v1,blue_sky_v1,blue_error_v1,transit_times_v1,exptime_v1],
                      names=['Net_photons','Gross_photons','Sky_photons','Error_photons',
                             'Time_hr','Exptime'])
        dat_r = Table([red_net_v1,red_gross_v1,red_sky_v1,red_error_v1,transit_times_v1,exptime_v1],
                      names=['Net_photons','Gross_photons','Sky_photons','Error_photons',
                             'Time_hr','Exptime'])

        ascii.write(dat_b,'visit1_blue.txt',overwrite=True)
        ascii.write(dat_r,'visit1_red.txt',overwrite=True)

    else:
        exptime_v2.append(exptime)
        transit_times_v2.append(transit_time)
        blue_error_v2.append(blue_err)
        red_error_v2.append(red_err)
        red_net_v2.append(red_net)
        blue_net_v2.append(blue_net)
        blue_gross_v2.append(blue_gross)
        red_gross_v2.append(red_gross)
        red_sky_v2 = np.array(red_gross_v2)-np.array(red_net_v2)
        blue_sky_v2 = np.array(blue_gross_v2)-np.array(blue_net_v2)

        dat_b = Table([blue_net_v2,blue_gross_v2,blue_sky_v2,blue_error_v2,transit_times_v2,exptime_v2],
                      names=['Net_photons','Gross_photons','Sky_photons','Error_photons',
                             'Time_hr','Exptime'])
        dat_r = Table([red_net_v2,red_gross_v2,red_sky_v2,red_error_v2,transit_times_v2,exptime_v2],
                      names=['Net_photons','Gross_photons','Sky_photons','Error_photons',
                             'Time_hr','Exptime'])

        ascii.write(dat_b,'visit2_blue.txt',overwrite=True)
        ascii.write(dat_r,'visit2_red.txt',overwrite=True)


'''
CALL WHICH PROCESS YOU WANT TO RUN
'''

GJ1132_TT = glob(dir_TT+'*TT1D*.fits')

exptime_v1 = []
exptime_v2 = []
blue_error_v1 = []
blue_error_v2 = []
red_error_v1 = []
red_error_v2 = []
transit_times_v1 = []
transit_times_v2 = []
red_net_v2 = []
red_net_v1 = []
blue_net_v2 = []
blue_net_v1 = []
blue_gross_v2 = []
blue_gross_v1 =[]
red_gross_v1 = []
red_gross_v2 = []

i = 0

for filename in GJ1132_TT:
    i += 1
    TT_FluxInt('GJ1132b', filename, 1215.8118, '10:14:51.77', '-47:09:24.1', 2457184.55786, 1.62893, 35.0)
