#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 07:59:41 2019

@author: irhamta
"""

import numpy as np
import pandas as pd
import os, fnmatch
import matplotlib.pyplot as plt
from scipy import interp


data = pd.read_csv('spex_list.csv')


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def response_func(template_spectrum, filter_response):

    '''
        This function can be used to calculate the AB-magnitude from spectral
        template multiplide with each bandpass response function.
    '''

    import numpy as np
    from scipy.integrate import simps

    c_AAs      = 2.99792458e18      # speed of light in Angstrom/s

    lamS, spec = template_spectrum  # Two columns with wavelength (in Angstrom) and flux density (in erg/s/cm2/AA)
    lamF, filt = filter_response    # Two columns with wavelength and response in the range [0,1]

    # this is correction for Vista filter
    # they gave transmission in scale of 0-100,
    # and wavelength should be timed by 10
    if np.max(filt) > 20:
        lamF*=10
        filt/=100

    # outside interpolation range, the value will be zero
    filt_int   = np.interp(lamS, lamF, filt, left=0, right=0) # Interpolate to common wavelength axis

    # if template is outside bandpass effective wavelength,
    # then assign 0 flux as non-detection
    if np.sum(filt_int) == 0:
        fnu = 0
        mAB = -np.inf
    elif np.sum(filt_int) > 0:
        I1         = simps(spec*filt_int*lamS, lamS)      # Denominator
        I2         = simps(filt_int/lamS, lamS)           # Numerator
        fnu        = I1/I2/c_AAs                          # Average flux density
        mAB        = -2.5*np.log10(fnu) - 48.6            # AB magnitude
    else:
        print ('Error! Something wrong....')

    return mAB, fnu

def calc_wise_mag(K, K_W1, W1_W2):

    '''
        Calculate WISE magnitude from colors table of Skrzypek et al. (2015).
        The W1 and W2 values are extrapolated from UKIDSS K-band.
        The input of K should be in AB-magnitude.
    '''

    # AB to Vega magnitude conversion for UKIDSS K-band
    K = K - 1.9

    W1 = K - K_W1
    W2 = W1 - W1_W2

    # Vega to AB magnitude conversion
    W1 = W1 + 2.699
    W2 = W2 + 3.339


    return W1, W2

def mag_to_flux(mag, cwave):

    '''
        Convert AB-magnitude to flux density as function of wavelength.
    '''

    fnu = 10**((mag + 48.6)/(-2.5)) # erg/s/Hz

    c_AAs   = 2.99792458e18         # speed of light in Angstrom/s
    flambda = fnu*c_AAs/(cwave**2)

    return flambda



# list of bandpasses
K_band = 'UKIRT_UKIDSS.K.dat'

# center wavelength of K, W1, W2, and W2 reddest part
cwave = [2.2060e+04, 3.3682e+04, 4.6179e+04, 5.36500e+04]

# the table from Skrzypek et al. (2015)
constant = pd.read_csv('mltdwarf_color.csv')
constant.index = constant['SpT']

index = 0
max_lamb = []

err_count = 0

for j in range(len(data)):
    name = data.loc[j]['Object Name'].replace(' ', '')

    path = find('*'+name+'*', '/home/irhamta/eazy-photoz/templates/SpeX/templates')


    try:
        mltdwarf = pd.read_csv(path[0], delimiter='\t', comment='#', header=None,
                               names=['wl', 'fl', 'fl_err'])

        mltdwarf['wl'] = mltdwarf*1e4

    except:
        print ('Template not found:', name)
        template = None
        err_count+=1
        continue


#    plt.figure()
#    plt.plot(mltdwarf['wl'], mltdwarf['fl'])
#    plt.title(data['SpT'].loc[j])
#    plt.xlabel(r'Wavelength [$\rm \AA$]')
#    plt.ylabel('Normalized Flux')
#    plt.savefig('figures/'+name+'.png')
#    plt.title
#    plt.close('all')





    # calculate K-band magnitude and flux
    mag_K, flux_K = response_func((mltdwarf['wl'], mltdwarf['fl']),
                                  np.loadtxt('filters/' + K_band, unpack=True))

    # check reddest wavelength of each template
    max_lamb.append(np.max(mltdwarf['wl']))

    # convert flux density from Hz^-1 to Angstrom^-1
    flux_K = mag_to_flux(mag_K, cwave[0])






    SpT = data['SpT'].loc[j]


    # calculate WISE magnitudes
    mag_W1, mag_W2 = calc_wise_mag(mag_K,
                                   constant['K_W1'].loc[SpT],
                                   constant['W1_W2'].loc[SpT])

    # ...and the fluxes
    flux_W1 = mag_to_flux(mag_W1, cwave[1])
    flux_W2 = mag_to_flux(mag_W2, cwave[2])

    # plot the output
    plt.figure(index)
    plt.plot(mltdwarf['wl'], mltdwarf['fl'], 'b-')
    plt.plot(cwave[0], flux_K, 'go')

    plt.plot(cwave[1], flux_W1, 'go')
    plt.plot(cwave[2], flux_W2, 'go')
    plt.plot([4775.6, 6129.5], [0, 0], 'go') # PS1 g, r

    wavelength = np.concatenate([ #[4775.6, 6129.5],
                                   mltdwarf['wl'],
                                  [cwave[1], cwave[2], cwave[3]] ])
    flambda = np.concatenate([ #[0, 0],
                                mltdwarf['fl'],
                               [flux_W1, flux_W2, 0] ])



    wavelength_interp = np.arange(min(wavelength), max(wavelength)+1, 1)
    flambda_interp = interp(wavelength_interp, wavelength, flambda,
                        left = 0, right=np.nan)

    plt.plot(wavelength, flambda, 'k-.')
    plt.plot(wavelength_interp, flambda_interp, 'r-.')


    plt.title(data['SpT'].loc[j])
    plt.xlabel(r'Wavelength [$\rm \AA$]')
    plt.ylabel('Normalized Flux')
    plt.savefig('figures/'+name+'.png')
    plt.title
    plt.close('all')


    output = pd.DataFrame({'# Wavelength_(Angstrom)': wavelength_interp,
                           'Flux_(cgs)'             : flambda_interp})

    output.to_csv('extrapolated/' + SpT + '_' + name +'.csv',
                  index = False, sep='\t', header=None)

    index+=1