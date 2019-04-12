from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import radio_beam
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from deconv_sources import deconv_srcs
import regions
from functools import reduce
import fnmatch
from astropy.coordinates import Angle, SkyCoord

def get_ind(names): #returns indices where sources are detected in all bands
    data = fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')
    flux_inds = [np.where(np.isnan(data['ap_flux_'+n]) == False)[0] for n in names]
    #fit_ind =  np.where(data['good_fit_flag'] == True)[0]
    ind1 = reduce(np.intersect1d, flux_inds)
    #ind = np.intersect1d(ind1, fit_ind)
    ind = ind1
    return ind


def plot_SEDs(names):
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7':339.7672758867*u.GHz}
    #freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7_hr':339.7672758867*u.GHz}
    data = fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')
    ind = get_ind(names)

    alpha = [1.5,2,2.5]
    freq_x = np.array([freqs[n].value for n in names])
    for i in ind:
        #fluxes = [data['gauss_amp_'+n][i] for n in names]
        #flux_err = [data['gauss_amp_err_'+n][i] for n in names]

        fluxes = [data['ap_flux_'+n][i] for n in names]
        flux_err = [data['ap_flux_err_'+n][i] for n in names]

        F2 = fluxes[1]
        nu2 = freq_x[1]
        if data['D_ID'][i] == 13:
            print(fluxes)
        
        plt.figure()
        plt.errorbar(freq_x, fluxes, yerr=flux_err, linestyle='', marker='o')
               
        for j in range(len(alpha)):
            F1 = F2*(freq_x/nu2)**alpha[j]
            plt.plot(freq_x, F1, linestyle='-', label=r'$\alpha = $'+str(alpha[j]))

        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('aperture flux (Jy)')
        #plt.ylabel('gaussian amplitude (Jy)')
        plt.xlabel('frequency (GHz)')
        plt.legend()
        plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/SEDs/SED_'+str(data['D_ID'][i])+'_apflux_bgfitted_apflux_final.png', dpi=300)

plot_SEDs(['B3', 'B6', 'B7'])
