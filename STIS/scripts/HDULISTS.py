#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:14:54 2017

@author: willwaalkes
"""

from glob import glob
import astropy.io.fits as f
import astropy.io.fits.hdu.hdulist as h
import os
import sys


'''
HDULISTS.py was developed for the purpose of dealing with TIME-TAG STIS data. After _tag files have been processed
through the TIME TAG analysis pipeline (inttag & calstis), _flt.fits and _x1d.fits files are created with extensions
equal to 3x the number of sub-exposures (here, that number is 10). 

Once these _flt.fits files have been produced, this program will separate the n-extension files into separate files,
1 file for each sub-exposure.
'''

dir_in = '/Users/willwaalkes/Desktop/Colorado/ColoradoResearch/GJ1132b/TIME-TAG/Calibration_Files/'
dir_out = '/Users/willwaalkes/Desktop/Colorado/ColoradoResearch/GJ1132b/TIME-TAG/Calibration_Files/Separated_Exposures/'

FLT_List = glob(dir_in+'*TT_flt.fits')
print FLT_List
for filename in FLT_List:
    
    hdu = f.open(filename)
    #data = hdu[1].data
    #header = hdu[1].header
    #print np.shape(data)
    #os.sys.exit()
    #w = data['wavelength'].flatten()
    #os.sys.exit()
    name = filename.replace('/Users/willwaalkes/Desktop/Colorado/ColoradoResearch/GJ1132b/TIME-TAG/Calibration_Files/', '')
    obs = name.replace('flt.fits', '')

    x1 = h.HDUList(hdus=[hdu[0],hdu[1],hdu[2],hdu[3]])
    x1.info()
    #x1.writeto(dir_out+obs+'Ex01_flt.fits',overwrite=True)

    x2 = h.HDUList(hdus=[hdu[0],hdu[4],hdu[5],hdu[6]])
    x2.info()
    #x2.writeto(dir_out+obs+'Ex02_flt.fits',overwrite=True)
    
    x3 = h.HDUList(hdus=[hdu[0],hdu[7],hdu[8],hdu[9]])
    x3.info()
    #x3.writeto(dir_out+obs+'Ex03_flt.fits',overwrite=True)
    
    x4 = h.HDUList(hdus=[hdu[0],hdu[10],hdu[11],hdu[12]])
    x4.info()
    #x4.writeto(dir_out+obs+'Ex04_flt.fits',overwrite=True)
    
    '''
    x5 = h.HDUList(hdus=[hdu[0],hdu[13],hdu[14],hdu[15]])
    x5.info()
    x5.writeto(dir_out+obs+'Ex05_flt.fits',overwrite=True)
    
    x6 = h.HDUList(hdus=[hdu[0],hdu[16],hdu[17],hdu[18]])
    x6.info()
    x6.writeto(dir_out+obs+'Ex06_flt.fits',overwrite=True)
    
    x7 = h.HDUList(hdus=[hdu[0],hdu[19],hdu[20],hdu[21]])
    x7.info()
    x7.writeto(dir_out+obs+'Ex07_flt.fits',overwrite=True)
    
    x8 = h.HDUList(hdus=[hdu[0],hdu[22],hdu[23],hdu[24]])
    x8.info()
    x8.writeto(dir_out+obs+'Ex08_flt.fits',overwrite=True)
    
    x9 = h.HDUList(hdus=[hdu[0],hdu[25],hdu[26],hdu[27]])
    x9.info()
    x9.writeto(dir_out+obs+'Ex09_flt.fits',overwrite=True)
    
    x10 = h.HDUList(hdus=[hdu[0],hdu[28],hdu[29],hdu[30]])
    x10.info()
    x10.writeto(dir_out+obs+'Ex10_flt.fits',overwrite=True)
    '''