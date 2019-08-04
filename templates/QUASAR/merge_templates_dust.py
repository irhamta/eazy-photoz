#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 20:17:31 2018

@author: irhamta
"""

# import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import interp
from astropy.table import Table
import os

# =============================================================================
# Open quasar templates
# =============================================================================

# Harris et al. (2016)
# BOSS quasars, 2.1 < z < 3.5
harris = Table.read('original/harris15_composite.fits')
# cut unnecessary wavelength
harris = harris[(harris['WAVE'] >= 912) & (harris['WAVE'] <= 3000)]


# Selsing et al. (2016)
# X-Shooter quasars, 1 < z < 2
selsing = np.loadtxt('original/Selsing2015.dat', skiprows=1)


# Jensen et al. (2016)
jensen = Table.read('original/Jensen2016.fits')
jensen = jensen[(jensen['BASEWAVE'] >= 912) & (jensen['BASEWAVE'] <= 3000)]

# =============================================================================
# Define spectral grid
# =============================================================================

# template from Chiara (priv. comm.)
banados_all = np.loadtxt('original/ps1sample_composite_may2016_extrapolated.spc')
banados_hlya = np.loadtxt('original/ps1sample_composite_may2016_high_lya_extrapolated.spc')
banados_llya = np.loadtxt('original/ps1sample_composite_may2016_low_lya_extrapolated.spc')

# use wavelength grid from Banados+16
wl = banados_all[:, 0]

# =============================================================================
# Extrapolation
# =============================================================================

# extrapolate Selsing fluxes to Banados wavelength
selsing_fl = interp(wl, selsing[:, 0], selsing[:, 1],
                    left = 0, right=np.nan)


# =============================================================================
# Scaling for Harris' data
# =============================================================================

# continuum definition from Jensen et al. (2016)
# http://iopscience.iop.org/article/10.3847/1538-4357/833/2/199/pdf
cont         = (wl >= 1440) & (wl <= 1480)
harris_fl    = interp(wl, harris['WAVE'], harris['FLUX'],
                   left = 0, right=np.nan)
harris_scale = np.nanmedian(selsing_fl[cont]/harris_fl[cont])

# rescale the flux to match Selsing's template
harris_fl[wl > 3000] = selsing_fl[wl > 3000]/harris_scale


# =============================================================================
# Scaling for Jensen's data
# =============================================================================

#cont2        = (wl >= 2160) & (wl <= 2230)
#jensen_fl    = interp(wl, jensen['BASEWAVE'], jensen['COMP11'],
#                   left = 0, right=np.nan)
#jensen_scale = np.nanmedian(selsing_fl[cont]/jensen_fl[cont])
#jensen_scale2 = np.nanmedian(selsing_fl[cont2]/jensen_fl[cont2])
#
## rescale the flux to match Selsing's template
#jensen_fl[wl > 2700] = selsing_fl[wl > 2700]/jensen_scale2


def stitch_spec(wl, fl, wl_grid=wl, selsing_fl=selsing_fl):
    '''
        This function combine any quasar template with Selsing's template.
        The fluxes will be stitched, rescaled at 2 continua.
    '''

    cont1 = (wl_grid >= 1440) & (wl_grid <= 1480)
    cont2 = (wl_grid >= 2160) & (wl_grid <= 2230)

    fl_int = interp(wl_grid, wl, fl, left = 0, right=np.nan)

    scale1 = np.nanmedian(selsing_fl[cont1]/fl_int[cont1])
    scale2 = np.nanmedian(selsing_fl[cont2]/fl_int[cont2])

    fl_int[wl_grid > 2650] = selsing_fl[wl_grid > 2650]/scale2

    return fl_int*scale1



plt.plot(wl, selsing_fl, 'b-', linewidth=2)
plt.plot(wl, harris_fl*harris_scale, 'r--', linewidth=2)


# =============================================================================
# Save templates grid
# =============================================================================

# save output grid
out_array_selsing       = np.array(banados_all)
out_array_selsing[:, 0] = wl
out_array_selsing[:, 1] = selsing_fl

out_array_harris        = np.array(banados_all)
out_array_harris[:, 0]  = wl
out_array_harris[:, 1]  = harris_fl*harris_scale


out_array_jensen        = np.array(banados_all)
out_array_jensen[:, 0]  = wl



for i in range(27):

    out_array_jensen[:, 1]  = stitch_spec(jensen['BASEWAVE'],
                                          jensen['COMP%i' %(i+1)])

    np.savetxt('nored/Jensen2016_extrapolated_%i.dat' %(i+1), out_array_jensen)
    plt.plot(out_array_jensen[:, 0], out_array_jensen[:, 1], '--')

plt.plot(wl, selsing_fl, 'k-')

np.savetxt('nored/Selsing2015_extrapolated.dat', out_array_selsing)
np.savetxt('nored/Harris2015_extrapolated.dat', out_array_harris)

plt.close('all')



# =============================================================================
# Add dust quasar spectral template
# =============================================================================

#ebv = np.array([0.00, 0.025, 0.05, 0.10])
ebv = np.arange(0, 0.142, 0.02)
Av = 4.05*ebv


def redden (wavelength, flux, Av):

    from extinction import fm07, apply, calzetti00

    # Fitzpatrick & Massa (2007) extinction model for R_V = 3.1
#    flux_red = apply(fm07(wavelength, Av), flux)

    # Calzetti (2000) extinction function.
    flux_red = apply(calzetti00(wave=wavelength, a_v=Av, r_v=4.05), flux)

    return flux_red


import pandas as pd
selsing_red = pd.DataFrame(out_array_selsing[:, 0], columns=['wave'])
harris_red = pd.DataFrame(out_array_harris[:, 0], columns=['wave'])

banados_all_red = pd.DataFrame(banados_all[:, 0], columns=['wave'])
banados_hlya_red = pd.DataFrame(banados_hlya[:, 0], columns=['wave'])
banados_llya_red = pd.DataFrame(banados_llya[:, 0], columns=['wave'])


tmp_fl = []

for av in Av:
    selsing_red[av/4.05] = redden(out_array_selsing[:, 0], out_array_selsing[:, 1], av)
    harris_red[av/4.05] = redden(out_array_harris[:, 0], out_array_harris[:, 1], av)

    banados_all_red[av/4.05] = redden(banados_all[:, 0], banados_all[:, 1], av)
    banados_hlya_red[av/4.05] = redden(banados_hlya[:, 0], banados_hlya[:, 1], av)
    banados_llya_red[av/4.05] = redden(banados_llya[:, 0], banados_llya[:, 1], av)

selsing_red.replace(np.nan, 0, inplace=True)
harris_red.replace(np.nan, 0, inplace=True)
banados_all_red.replace(np.nan, 0, inplace=True)
banados_hlya_red.replace(np.nan, 0, inplace=True)
banados_llya_red.replace(np.nan, 0, inplace=True)

plt.figure()
for col in selsing_red.columns[1:]:
    plt.plot(selsing_red['wave'], selsing_red[col])
#    plt.plot(harris_red['wave'], harris_red[col])

#    plt.plot(banados_all_red['wave'], banados_all_red[col])
#    plt.plot(banados_hlya_red['wave'], banados_hlya_red[col])
#    plt.plot(banados_llya_red['wave'], banados_llya_red[col])


    selsing_red[['wave', col]].to_csv('dusty/Selsing2015_extrapolated_calz00_ebv%.3f.dat' %col,
                                       header=False, index=False, sep='\t')
    harris_red[['wave', col]].to_csv('dusty/Harris2015_extrapolated_calz00_ebv%.3f.dat' %col,
                                       header=False, index=False, sep='\t')

    banados_all_red[['wave', col]].to_csv('dusty/Banados2016_avg_extrapolated_calz00_ebv%.3f.dat' %col,
                                       header=False, index=False, sep='\t')
    banados_hlya_red[['wave', col]].to_csv('dusty/Banados2016_hlya_extrapolated_calz00_ebv%.3f.dat' %col,
                                       header=False, index=False, sep='\t')
    banados_llya_red[['wave', col]].to_csv('dusty/Banados2016_llya_extrapolated_calz00_ebv%.3f.dat' %col,
                                       header=False, index=False, sep='\t')
plt.xlabel(r'Wavelength [$\rm \AA$]')
plt.ylabel(r'Normalized Flux')
plt.close('all')



# =============================================================================
# Apply reddening for Jensen+16 templates
# =============================================================================

jensen_list = [x for x in os.listdir('nored/') if 'Jensen' in x]

for file in jensen_list:
    jensen_tmp = np.loadtxt('nored/' + file)

    for av in Av:
        jensen_tmp[:, 1] = redden(jensen_tmp[:, 0], jensen_tmp[:, 1], av)

        # replace nan values
        jensen_tmp[np.isnan(jensen_tmp)] = 0

        np.savetxt('dusty/'+file[:-4]+'_calz00_ebv%.3f.dat' %(av/4.05), jensen_tmp)

        plt.plot(jensen_tmp[:, 0], jensen_tmp[:, 1])

