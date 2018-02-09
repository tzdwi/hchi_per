"""
We already have done all the calibration steps. Now it's just plugging and chugging on the data, with some degree of interactivity. Let's go.

Author: Trevor Dorn-Wallenstein
Date: 2/7/18
"""

import pydis
from pydis.wrappers import _WriteSpec as WriteSpec
import numpy as np
from glob import glob
from astropy.io import fits

data_dir = '../data/1_25_18/Q1UW03/UT180126/'

#Load in our calibration results
#biases
bbias = fits.getdata(data_dir+'BIAS_B.fits')
rbias = fits.getdata(data_dir+'BIAS_R.fits')

#Flats and masks for where the response is good
bflat = fits.getdata(data_dir+'FLAT_B.fits')
bfmask_out = np.genfromtxt(data_dir+'FMASK_B.txt')
rflat = fits.getdata(data_dir+'FLAT_R.fits')
rfmask_out = np.genfromtxt(data_dir+'FMASK_R.txt')

#Wavelength solution from the night
bwfit = np.genfromtxt(data_dir+'WAVE_BLUE.txt')
rwfit = np.genfromtxt(data_dir+'WAVE_RED.txt')

#Sensititivty function from the night
bsensfunc = np.genfromtxt(data_dir+'SENS_B.txt')
rsensfunc = np.genfromtxt(data_dir+'SENS_R.txt')

print('Loaded in calibration data')
