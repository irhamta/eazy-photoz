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
banados_all = np.loadtxt('ps1sample_composite_may2016_extrapolated.spc')

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


for i in [18]:#range(27):

    out_array_jensen[:, 1]  = stitch_spec(jensen['BASEWAVE'],
                                          jensen['COMP%i' %(i+1)])

    np.savetxt('Jensen2016_extrapolated_%i.dat' %(i+1), out_array_jensen)
    plt.plot(out_array_jensen[:, 0], out_array_jensen[:, 1], '--')

#plt.plot(wl, selsing_fl, 'k-')

np.savetxt('Selsing2015_extrapolated.dat', out_array_selsing)
np.savetxt('Harris2015_extrapolated.dat', out_array_harris)

plt.close('all')


#test = Table.read('../MBrown_ATLAS/NGC_6240_spec.dat', format='ascii')
#plt.plot(test['col1'], test['col2'])
#brak

banados_hlya = np.loadtxt('ps1sample_composite_may2016_high_lya_extrapolated.spc')
banados_llya = np.loadtxt('ps1sample_composite_may2016_low_lya_extrapolated.spc')





# continuum to scale templates
cont = (wl >= 1285) & (wl <= 1295)


banados_all_scale = np.nanmedian(selsing_fl[cont]/banados_all[:, 1][cont])
banados_hlya_scale = np.nanmedian(selsing_fl[cont]/banados_hlya[:, 1][cont])
banados_llya_scale = np.nanmedian(selsing_fl[cont]/banados_llya[:, 1][cont])

plt.plot(wl, selsing_fl, 'b-')
plt.plot(wl, banados_all[:, 1]*banados_all_scale, 'r-')
plt.plot(wl, banados_hlya[:, 1]*banados_hlya_scale, 'g-')
plt.plot(wl, banados_llya[:, 1]*banados_llya_scale, 'y-')
#plt.close('all')



#cont2 = (wl >= 1340) & (wl <= 1360)
#jensen_fl = interp(wl, jensen['BASEWAVE'], jensen['COMP18'], left = 0, right=np.nan)
#jensen_scale = np.nanmedian(selsing_fl[cont]/jensen_fl[cont])
#jensen_fl[wl > 3000] = selsing_fl[wl > 3000]/jensen_scale
#plt.plot(wl, jensen_fl*jensen_scale, 'r--')




#brak
#for j in range(27):
#    plt.plot(jensen['BASEWAVE'], jensen['COMP%i' %(j+1)])

plt.xlim(0, 3000)

#plt.plot(banados[:, 0], banados_hlya[:, 1]*median_flux, 'g-')
#plt.plot(banados[:, 0], banados_llya[:, 1]*median_flux, 'y-')