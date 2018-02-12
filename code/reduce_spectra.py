"""
We already have done all the calibration steps. Now it's just plugging and chugging on the data, with some degree of interactivity. Let's go.

Author: Trevor Dorn-Wallenstein
Date: 2/7/18
"""

import numpy as np
import pydis
from pydis.wrappers import _WriteSpec as WriteSpec
import os
from glob import glob
from astropy.io import fits
from os import remove
from matplotlib import pyplot as plt, cm

data_dir = '/Volumes/shoobert/Research/UW/hchi_per/data/1_25_18/Q1UW03/UT180126/'

#ask if we want to do calibrations

#biases
biascheck = ''
while biascheck not in ['y','n']:
    biascheck = input("Use old biases? [y/n]: ")
    
if biascheck == 'n':
    bbias = pydis.biascombine(data_dir+'bias_b.lis', trim=True, output = data_dir+'BIAS_B.fits')
    rbias = pydis.biascombine(data_dir+'bias_r.lis', trim=True, output = data_dir+'BIAS_R.fits')
    
elif biascheck == 'y':
    bbias = fits.getdata(data_dir+'BIAS_B.fits')
    rbias = fits.getdata(data_dir+'BIAS_R.fits')
    
#flats
flatcheck = ''
while flatcheck not in ['y','n']:
    flatcheck = input("Use old flats? [y/n]: ")
    
if flatcheck == 'n':
    mode = input("spline or polynomial? [s/p]: ")
    order = int(input("order of spline/poly? [int]: "))
    if mode == 's':
        bflat,bfmask_out = pydis.flatcombine(data_dir+'flat_b.lis', bbias, mode='spline', spline_order=order, output=data_dir+'FLAT_B.fits')
        rflat,rfmask_out = pydis.flatcombine(data_dir+'flat_r.lis', rbias, mode='spline', spline_order=order, output=data_dir+'FLAT_R.fits')
    elif mode == 'p':
        bflat,bfmask_out = pydis.flatcombine(data_dir+'flat_b.lis', bbias, flat_poly=order, output=data_dir+'FLAT_B.fits')
        rflat,rfmask_out = pydis.flatcombine(data_dir+'flat_r.lis', rbias, flat_poly=order, output=data_dir+'FLAT_R.fits')
     
    np.savetxt(data_dir+'FMASK_B.txt',bfmask_out)
    np.savetxt(data_dir+'FMASK_R.txt',rfmask_out)
    
elif flatcheck == 'y':
    bflat = fits.getdata(data_dir+'FLAT_B.fits')
    bfmask_out = np.genfromtxt(data_dir+'FMASK_B.txt').astype(int)
    rflat = fits.getdata(data_dir+'FLAT_R.fits')
    rfmask_out = np.genfromtxt(data_dir+'FMASK_R.txt').astype(int)
    
    
#wavelength calibration
wavecheck = ''
while wavecheck not in ['y','n']:
    wavecheck = input("Use old wavelength solution? [y/n]: ")
    
if wavecheck == 'n':
    bwfit = pydis.HeNeAr_fit(data_dir+'HeNeAr.0005b.fits', trim=True, fmask=bfmask_out,interac=True, mode='poly', fit_order=5)
    rwfit = pydis.HeNeAr_fit(data_dir+'HeNeAr.0005r.fits', trim=True, fmask=rfmask_out,interac=True, mode='poly', fit_order=5)
    
    np.savetxt('WAVE_BLUE.txt',bwfit)
    np.savetxt('WAVE_RED.txt',rwfit)
    
elif wavecheck == 'y':
    bwfit = np.genfromtxt(data_dir+'WAVE_BLUE.txt')
    rwfit = np.genfromtxt(data_dir+'WAVE_RED.txt')
    
senscheck = ''
while senscheck not in ['y','n']:
    senscheck = input("Use old sensitivity function? [y/n]: ")
    
if senscheck == 'n':
    
    bstdlist = glob(data_dir+'bd52*b.fits')
    rstdlist = glob(data_dir+'bd52*r.fits')
    
    for stdlist,col in zip([bstdlist,rstdlist],['b','r']):
        if col == 'b':
            bias = bbias
            flat = bflat
            fmask_out = bfmask_out
            wfit = bwfit
        elif col == 'r':
            bias = rbias
            flat = rflat
            fmask_out = rfmask_out
            wfit = rwfit
        
        #trace the first one so we have a wavelength grid to interpolate to, and define our apertures
        stdspec = stdlist[0]
        std = pydis.OpenImg(stdspec, trim=True)

        stdraw = std.data
        stdexptime = std.exptime
        stdairmass = std.airmass
        stddata = ((stdraw - bias) / flat) / stdexptime
        
        stdtrace = pydis.ap_trace(stddata, fmask=fmask_out, nsteps=15, interac=False, display=True)
        wfinalstd = pydis.mapwavelength(stdtrace, wfit, mode='poly')
        
        print('Defining apertures...')
        done = 'n'
        xbins = np.arange(stddata.shape[1])
        while done == 'n':
            
            apwidth = int(input('apwidth = '))
            skysep = int(input('skysep = '))
            skywidth = int(input('skywidth = '))

            plt.figure()
            plt.imshow(np.log10(stddata), origin='lower',aspect='auto',cmap=cm.Greys_r)

            # the trace
            plt.plot(xbins, stdtrace,'b',lw=1)

            # the aperture
            plt.plot(xbins, stdtrace-apwidth,'r',lw=1)
            plt.plot(xbins, stdtrace+apwidth,'r',lw=1)

            # the sky regions
            plt.plot(xbins, stdtrace-apwidth-skysep,'g',lw=1)
            plt.plot(xbins, stdtrace-apwidth-skysep-skywidth,'g',lw=1)
            plt.plot(xbins, stdtrace+apwidth+skysep,'g',lw=1)
            plt.plot(xbins, stdtrace+apwidth+skysep+skywidth,'g',lw=1)

            plt.title('(with trace, aperture, and sky regions)')
            plt.show()
            
            done = input('Done? [y/n]: ')
        
        stdfluxes = []
        for stdspec in stdlist:
            stdimg = pydis.OpenImg(stdspec, trim=True)
            stdraw = stdimg.data
            stdexptime = stdimg.exptime
            stdairmass = stdimg.airmass
            stddata = ((stdraw - bias) / flat) / stdexptime
            stdtrace = pydis.ap_trace(stddata, fmask=fmask_out, nsteps=15, interac=False, display=True)
            ext_std, stdsky, stderr = pydis.ap_extract(stddata, stdtrace, apwidth=apwidth, skysep=skysep, skywidth=skywidth, skydeg=0)
            std_flux_temp = ext_std - stdsky
            wfinalstd_tmp = pydis.mapwavelength(stdtrace, wfit, mode='poly')
            std_flux_temp_x = pydis.AirmassCor(wfinalstd_tmp,std_flux_temp,stdairmass, airmass_file='apoextinct.dat')
            if col == 'b':
                std_flux = np.interp(wfinalstd,wfinalstd_tmp,std_flux_temp_x)
            elif col == 'r':
                std_flux = np.interp(wfinalstd[::-1],wfinalstd_tmp[::-1],std_flux_temp_x[::-1])
                std_flux = std_flux[::-1]
                
            stdfluxes.append(std_flux)
            
        stdflux = np.median(stdfluxes,axis=0)
        sensfunc = pydis.DefFluxCal(wfinalstd, stdflux, mode='spline',stdstar='spec50cal/g191b2b.dat')
        if col == 'b':
            bsensfunc = sensfunc
            np.savetxt(data_dir+'SENS_B.txt',sensfunc)
            np.savetxt(data_dir+'WFINALSTD_B.txt',wfinalstd)
        if col == 'r':
            rsensfunc = sensfunc
            np.savetxt(data_dir+'SENS_R.txt',sensfunc)
            np.savetxt(data_dir+'WFINALSTD_R.txt',wfinalstd)
        
elif senscheck == 'y':
    bsensfunc = np.genfromtxt(data_dir+'SENS_B.txt')
    rsensfunc = np.genfromtxt(data_dir+'SENS_R.txt')
    bwfinalstd = np.genfromtxt(data_dir+'WFINALSTD_B.txt')
    rwfinalstd = np.genfromtxt(data_dir+'WFINALSTD_B.txt')
    

print('Finished with calibration data')

do_reductions = input("Are you ready to slog through the data? [y/n]: ")

if do_reductions == 'y':

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
                wfinalstd = bwfinalstd
            elif col == 'r':
                bias = rbias
                flat = rflat
                fmask_out = rfmask_out
                wfit = rwfit
                sensfunc = rsensfunc
                wfinalstd = rwfinalstd
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

            plt.plot(wfinal,flux)
            plt.show()