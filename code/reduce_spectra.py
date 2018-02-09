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
from os import remove

data_dir = '../data/1_25_18/Q1UW03/UT180126/'

#Load in our calibration results
#biases
bbias = fits.getdata(data_dir+'BIAS_B.fits')
rbias = fits.getdata(data_dir+'BIAS_R.fits')

#Flats and masks for where the response is good
bflat = fits.getdata(data_dir+'FLAT_B.fits')
bfmask_out = np.genfromtxt(data_dir+'FMASK_B.txt').astype(int)
rflat = fits.getdata(data_dir+'FLAT_R.fits')
rfmask_out = np.genfromtxt(data_dir+'FMASK_R.txt').astype(int)

#Wavelength solution from the night
bwfit = np.genfromtxt(data_dir+'WAVE_BLUE.txt')
rwfit = np.genfromtxt(data_dir+'WAVE_RED.txt')

#Sensititivty function from the night
bsensfunc = np.genfromtxt(data_dir+'SENS_B.txt')
rsensfunc = np.genfromtxt(data_dir+'SENS_R.txt')

#wavelength array from sensfunc calculation
bwfitstd = np.genfromtxt(data_dir+'WFITSTD_B.txt')
rwfitstd = np.genfromtxt(data_dir+'WFITSTD_B.txt')

print('Loaded calibration data')

fitsfiles = glob(data_dir+'*.fits')

for file in fitsfiles:
    
    name = file.split('/')[-1]
    
    #Check if the file is a science target
    if name[:3] in ['blu','red','hei']:
        print('Reducing {}'.format(name))
        #Give me the appropriate cal data based on if its a red or blue exposure
        col = name.split('.')[1][-1]
        if col == 'b':
            bias = bbias
            flat = bflat
            fmask_out = bfmask_out
            wfit = bwfit
            sensfunc = bsensfunc
            wfitstd = bwfitstd
        elif col == 'r':
            bias = rbias
            flat = rflat
            fmask_out = rfmask_out
            wfit = rwfit
            sensfunc = rsensfunc
            wfitstd = rwfitstd
        else:
            print('something went wrong with file parsing')
            break
            
        # open the image, get the data and a couple important numbers
        print('Reading in and reducing image...')
        img = pydis.OpenImg(file, trim=True)
        raw = img.data
        exptime = img.exptime
        airmass = img.airmass

        # remove bias and flat, divide by exptime
        data = ((raw - bias) / flat) / exptime
        
        # Do the trace
        print('Tracing spectrum...')
        done = 'n'
        while done == 'n':
            # trace the science image
            nsteps = int(input('nsteps = '))
            trace = pydis.ap_trace(data, fmask=fmask_out, nsteps=nsteps, interac=True, display=True)
            done = input('Done? [y/n]: ')
        
        # Now do the extraction region defining
        print('Defining apertures...')
        done = 'n'
        xbins = np.arange(data.shape[1])
        while done == 'n':
            
            apwidth = int(input('apwidth = '))
            skysep = int(input('skysep = '))
            skywidth = int(input('skywidth = '))

            plt.figure()
            plt.imshow(np.log10(data), origin='lower',aspect='auto',cmap=cm.Greys_r)

            # the trace
            plt.plot(xbins, trace,'b',lw=1)

            # the aperture
            plt.plot(xbins, trace-apwidth,'r',lw=1)
            plt.plot(xbins, trace+apwidth,'r',lw=1)

            # the sky regions
            plt.plot(xbins, trace-apwidth-skysep,'g',lw=1)
            plt.plot(xbins, trace-apwidth-skysep-skywidth,'g',lw=1)
            plt.plot(xbins, trace+apwidth+skysep,'g',lw=1)
            plt.plot(xbins, trace+apwidth+skysep+skywidth,'g',lw=1)

            plt.title('(with trace, aperture, and sky regions)')
            plt.show()
            
            done = input('Done? [y/n]: ')
            
        
        #Extraction and sky subtraction
        print('Extracting and sky subtracting...')
        ext_spec, sky, fluxerr = pydis.ap_extract(data, trace, apwidth=apwidth, skysep=skysep, skywidth=skywidth, skydeg=0)
        flux = ext_spec - sky
        
        #Calibrating and writing...
        print('Calibrating and writing to {0}.spec and {0}.raw.spec'.format(name))
        
        wfinal = pydis.mapwavelength(trace, wfit, mode='poly')
        flux_x = pydis.AirmassCor(wfinal, flux, airmass, airmass_file='apoextinct.dat')
        ffinal,efinal = pydis.ApplyFluxCal(wfinal, flux_x, fluxerr, wfinalstd, sensfunc)
        
        WriteSpec(file+'.raw', wfinal, flux, fluxerr, trace)
        os.remove(file+'.raw.trace')
        WriteSpec(file, wfinal, ffinal, efinal, trace)