#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:13:53 2017

@author: willwaalkes
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import batman
import scipy.stats as stats

blue_shift_v1 = np.loadtxt('blueshiftv1.txt')
blue_shift_v1 = np.array(blue_shift_v1)
red_shift_v1 = np.loadtxt('redshiftv1.txt')
red_shift_v1 = np.array(red_shift_v1)
blue_error_v1 = np.loadtxt('blueerrorv1.txt')
blue_error_v1 = np.array(blue_error_v1)
red_error_v1 = np.loadtxt('rederrorv1.txt')
red_error_v1 = np.array(red_error_v1)
transit_times_v1 = np.loadtxt('transittimesv1.txt')
transit_times_v1 = np.array(transit_times_v1)
blue_shift_v2 = np.loadtxt('blueshiftv2.txt')
blue_shift_v2 = np.array(blue_shift_v2)
red_shift_v2 = np.loadtxt('redshiftv2.txt')
red_shift_v2 = np.array(red_shift_v2)
blue_error_v2 = np.loadtxt('blueerrorv2.txt')
blue_error_v2 = np.array(blue_error_v2)
red_error_v2 = np.loadtxt('rederrorv2.txt')
red_error_v2 = np.array(red_error_v2)
transit_times_v2 = np.loadtxt('transittimesv2.txt')
transit_times_v2 = np.array(transit_times_v2)

desktop = '/Users/willwaalkes/Desktop/'

def BATMAN_MODEL(Rp, Baseline, t = None):
    #Next we create a TransitParams object to store the physical parameters describing the transit:

    params = batman.TransitParams()
    params.t0 = 0.                       #time of inferior conjunction
    params.per = 39.09432                #period in hours
    params.rp = Rp                       #planet radius (in units of stellar radii)
    params.a = 16.                       #semi-major axis (in units of stellar radii)
    params.inc = 90.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.1956, 0.3700]          #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model

    m = batman.TransitModel(params, t)    #initializes model

    flux = m.light_curve(params)*Baseline        #calculates light curve
    return flux

Rps = np.linspace(0,3,151)
Baselines = np.linspace(0.0,1.4,141)

highres_times = np.linspace(-18, 2.5, 500)

'''
write a function that takes in the rp and baseline, and spits out the chi2 value
feed in different values for the baseline, and the function will spit out all the chi2

Then plot chi2 on the y axis and baseline on the x axis

make a 2 D plot with a chi2 colorbar or contour

plt.scatter with colors assigned to points (x,y,c=z, cmap = )
'''

def fixed_rp(times, baseline, rp, data, error):

    chi = np.zeros(shape=(len(data),len(baseline)))

    '''
    first, run the BATMAN model

    second, calculate the total chi-square for the model, data
    '''
    for i in range(len(baseline)):
        flux = BATMAN_MODEL(rp, baseline[i], times)

        for j in range(len(flux)):
            chi2 = ((data[j]-flux[j])/error[j])**2
            chi[j,i] = chi2

    chi2sum = np.sum(chi, axis = 0)

    '''
    Last, plot the chi2 (x axis) vs baselines (y axis)
    '''

    plt.plot(chi2sum,baseline)
    plt.ylabel('Baseline')
    plt.xlabel('Chi Square')
    plt.show()
    plt.clf()

def fixed_base(times, baseline, rp, data, error):

    chi = np.zeros(shape=(len(data),len(rp)))

    '''
    First, run the BATMAN model

    Second, calculate the total chi-square for the model, data
    '''
    for i in range(len(rp)):
        flux = BATMAN_MODEL(rp[i], baseline, times)

        for j in range(len(flux)):
            chi2 = ((data[j]-flux[j])/error[j])**2
            chi[j,i] = chi2

    chi2sum = np.sum(chi, axis = 0)

    '''
    Last, plot the chi2 (x axis) vs rp (y axis)
    '''

    plt.plot(chi2sum,rp)
    plt.ylabel('Rp')
    plt.xlabel('Chi Square')
    plt.show()
    plt.clf()

def all_params(times, baseline, rp, data, error, mask, title):

    chi = np.zeros(shape=(len(data),len(rp),len(baseline)))

    '''
    Calculate the scatter
    '''

    times = np.delete(times,mask)
    data = np.delete(data,mask)
    error = np.delete(error,mask)

    Best_Base = np.sum(data/error**(2.0))/np.sum(1/error**(2.0))
    sigma_sq = np.sum((data-Best_Base)**2)/(len(data)-1.0)
    error_avg = np.sum(error**(2.0))/len(error)
    sigma_scat = np.sqrt(sigma_sq-error_avg)

    error_tot = np.sqrt(sigma_scat**2+error**2)

    for i in range(len(rp)):
        for j in range(len(baseline)):

            flux = BATMAN_MODEL(rp[i], baseline[j], times)

            for k in range(len(flux)):
                chi2 = ((data[k]-flux[k])/error_tot[k])**2
                chi[k,i,j] = chi2

    chi2sum = np.sum(chi, axis = 0)

    p_value = 1 - stats.chi2.cdf(x=chi2sum,df=28)

    chi2_min = np.amin(chi2sum)
    indices = np.where(chi2sum == chi2_min)
    print chi2_min/28
    print indices
    #os.sys.exit()
    Best_RP = rp[indices[0]]
    print "Best RP",Best_RP
    Best_Base = baseline[indices[1]]
    print "Best Base",Best_Base
    print '------'

    highres_model = BATMAN_MODEL(Best_RP, Best_Base, t=highres_times)


    plt.plot(highres_times, highres_model,color='r')#,linewidth=0.15*highres_model)
    plt.errorbar(highres_times, highres_model,yerr=0.15*highres_model,color='gray')
    plt.errorbar(times,data,yerr=error,fmt='o',color='b')
    plt.xlabel('Time (Hours from Mid-Transit)')
    plt.ylabel('Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$)')
    plt.title(Best_RP,loc='center')
    plt.title(Best_Base,loc='right')
    plt.title('Best Fit', loc='left')
    plt.show()    

    plt.imshow(chi2sum[:,:],extent=(baseline[0],baseline[-1],rp[-1],rp[0]),aspect='auto',cmap='coolwarm',
               vmin=chi2_min,vmax=(chi2_min+10.0))
    plt.xlabel('Baseline Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$)')
    plt.ylabel('Rp')
    plt.colorbar(label='$\chi^2$')
    plt.title(title)
    plt.show()
    plt.clf()

    '''
    Fixed_Base = chi2sum[:,indices[1]]
    plt.plot(Rps,Fixed_Base)
    plt.axhline(chi2_min,color='purple')
    plt.axhline(chi2_min+1,color='g')
    plt.axhline(chi2_min+4,color='r')
    plt.xlim(0.1,0.4)
    plt.ylim(chi2_min,chi2_min+4.5)
    plt.xlabel('Rps')
    plt.ylabel('Chi2')
    plt.title('Fixed Base')
    plt.show()
    plt.clf()

    Fixed_RP = chi2sum[indices[0],:]
    plt.plot(Baselines,Fixed_RP[0,:])
    plt.axhline(chi2_min,color='purple')
    plt.axhline(chi2_min+1,color='g')
    plt.axhline(chi2_min+4,color='r')
    plt.ylim(chi2_min,chi2_min+4.5)
    plt.xlim(1.1,1.3)
    plt.xlabel('Baselines')
    plt.ylabel('Chi2')
    plt.title('Fixed RP')
    plt.show()
    plt.clf()
    '''

    plt.imshow(np.log10(p_value[:,:]),extent=(baseline[0],baseline[-1],rp[-1],rp[0]),aspect='auto',
               cmap='coolwarm',vmin=-5)
    plt.xlabel('Baseline Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$)')
    plt.ylabel('Rp')
    plt.colorbar(label='log$_{10}$(p)')
    plt.title('p-values')
    plt.show()
    plt.clf()


def combined_plots(times, baseline, rp, data_blue, data_red, error_blue, sig_scat, error_red,title):

    chi_blue = np.zeros(shape=(len(data_blue),len(rp),len(baseline)))
    chi_red = np.zeros(shape=(len(data_red),len(rp),len(baseline)))

    '''
    Blue wing calculations
    '''

    Best_Base_blue = np.sum(data_blue/error_blue**(2.0))/np.sum(1/error_blue**(2.0))
    sigma_sq_blue = np.sum((data_blue-Best_Base_blue)**2)/(len(data_blue)-1.0)
    error_avg_blue = np.sum(error_blue**(2.0))/len(error_blue)
    sigma_scat_blue = np.sqrt(sigma_sq_blue-error_avg_blue)

    error_tot_blue = np.sqrt(sigma_scat_blue**2+error_blue**2)

    for i in range(len(rp)):
        for j in range(len(baseline)):

            flux = BATMAN_MODEL(rp[i], baseline[j], times)

            for k in range(len(flux)):
                chi2_blue = ((data_blue[k]-flux[k])/error_tot_blue[k])**2
                chi_blue[k,i,j] = chi2_blue

    chi2sum_blue = np.sum(chi_blue, axis = 0)

    chi2_min_blue = np.amin(chi2sum_blue)
    indices_blue = np.where(chi2sum_blue == chi2_min_blue)
    print chi2_min_blue
    print indices_blue

    Best_RP_blue = rp[indices_blue[0]]
    print "blue RP",Best_RP_blue
    Best_Base_blue = baseline[indices_blue[1]]
    print "blue Base",Best_Base_blue

    highres_model_blue = BATMAN_MODEL(Best_RP_blue, Best_Base_blue, t=highres_times)

    '''
    Red wing calculations
    '''

    Best_Base_red = np.sum(data_red/error_red**(2.0))/np.sum(1/error_red**(2.0))
    sigma_sq_red = np.sum((data_red-Best_Base_red)**2)/(len(data_red)-1.0)
    error_avg_red = np.sum(error_red**(2.0))/len(error_red)
    sigma_scat_red = np.sqrt(sigma_sq_red-error_avg_red)

    error_tot_red = np.sqrt(sigma_scat_red**2+error_red**2)

    for i in range(len(rp)):
        for j in range(len(baseline)):

            flux = BATMAN_MODEL(rp[i], baseline[j], times)

            for k in range(len(flux)):
                chi2_red = ((data_red[k]-flux[k])/error_tot_red[k])**2
                chi_red[k,i,j] = chi2_red

    chi2sum_red = np.sum(chi_red, axis = 0)

    chi2_min_red = np.amin(chi2sum_red)
    indices_red = np.where(chi2sum_red == chi2_min_red)
    print chi2_min_red
    print indices_red

    Best_RP_red = rp[indices_red[0]]
    print "Red RP",Best_RP_red
    Best_Base_red = baseline[indices_red[1]]
    print "Red Base",Best_Base_red

    highres_model_red = BATMAN_MODEL(Best_RP_red, Best_Base_red, t=highres_times)

    '''
    PLOT BOTH ON THE SAME AXIS
    '''

    plt.plot(highres_times, highres_model_blue,color='royalblue',zorder=10)
    plt.plot(highres_times, highres_model_red,color='crimson',zorder=10)
    plt.errorbar(highres_times, highres_model_blue,yerr=sig_scat*highres_model_blue,color='gray',alpha=0.5,zorder=-100)
    plt.errorbar(highres_times, highres_model_red,yerr=sig_scat*highres_model_red,color='gray',alpha=0.5,zorder=-100)
    plt.errorbar(times,data_blue,yerr=error_blue,fmt='o',color='royalblue',zorder=100)
    plt.errorbar(times,data_red,yerr=error_red,fmt='o',color='crimson',zorder=100)
    plt.xlabel('Time (Hours from Mid-Transit)',fontsize=15)
    plt.ylabel('Flux ($x10^{-14}$ erg s$^{-1}$ cm$^{-2}$)',fontsize=15)
    plt.ylim(-0.4,1.7)
    plt.title(title,fontsize=18)
    plt.tight_layout()
    plt.savefig(desktop+title+'_lightcurve.pdf',dpi=500)
    plt.show()
    plt.clf()

def variability(data, error):

    Best_Base = np.sum(data/error**(2.0))/np.sum(1/error**(2.0))

    sigma_sq = np.sum((data-Best_Base)**2)/(len(data)-1.0) #This gives a single value for sigma

    error_avg = np.sum(error**(2.0))/len(error)

    sigma_scat = np.sqrt(sigma_sq-error_avg)

    variability = sigma_scat/Best_Base
    print 'variability',variability
    print "min variability:",np.nanmin(variability)
    print "max variability:",np.nanmax(variability)
    print ''


#all_params(times=transit_times_v1, baseline=Baselines, rp=Rps, data=blue_shift_v1, error=blue_error_v1, mask = [0,6,10,22], title='$\chi^2$ Map for Visit 1 Blue Wing')
#all_params(times=transit_times_v1, baseline=Baselines, rp=Rps, data=red_shift_v1, error=red_error_v1, mask = [0,6,10,22], title='$\chi^2$ Map for Visit 1 Red Wing')
#all_params(times=transit_times_v2, baseline=Baselines, rp=Rps, data=blue_shift_v2, error=blue_error_v2, mask = [9,13,18,21], title='$\chi^2$ Map for Visit 2 Blue Wing')
#all_params(times=transit_times_v2, baseline=Baselines, rp=Rps, data=red_shift_v2, error=red_error_v2, mask = [9,13,18,21], title='$\chi^2$ Map for Visit 2 Red Wing')

combined_plots(times=transit_times_v1, baseline=Baselines, rp=Rps,
               data_blue=blue_shift_v1, data_red=red_shift_v1,
               error_blue=blue_error_v1, error_red=red_error_v1,
               sig_scat = 0.14, title='Visit 1')
combined_plots(times=transit_times_v2, baseline=Baselines, rp=Rps,
               data_blue=blue_shift_v2, data_red=red_shift_v2,
               error_blue=blue_error_v2, error_red=red_error_v2,
               sig_scat = 0.09, title='Visit 2')

#variability(data=blue_shift_v1, error=blue_error_v1)
#variability(data=red_shift_v1, error=red_error_v1)
#variability(data=blue_shift_v2, error=blue_error_v2)
#variability(data=red_shift_v2, error=red_error_v2)
