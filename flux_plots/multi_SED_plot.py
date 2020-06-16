from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from functools import reduce
from astropy.coordinates import Angle, SkyCoord

def get_ind(table): #returns indices where sources are detected in all bands
    names = ['B3', 'B6', 'B7']
    flux_inds = [np.where(np.isnan(table['ap_flux_'+n]) == False)[0] for n in names]
    ind = reduce(np.intersect1d, flux_inds)
    
    return ind


def plot_SED(ax, srctab):
    freqs = np.array([98, 223.5, 339.7672758867]) #*u.GHz
    
    alpha = [1.5,2,2.5]

    fluxes = [srctab['ap_flux_B3'], srctab['ap_flux_B6'], srctab['ap_flux_B7']]
    flux_errs = [srctab['ap_flux_err_B3'], srctab['ap_flux_err_B6'], srctab['ap_flux_err_B7']]

    F2 = fluxes[1]
    nu2 = freqs[1]
   
    ax.errorbar(freqs, fluxes, yerr=flux_errs, linestyle='', marker='o')
               
    for j in range(len(alpha)):
        F1 = F2*(freqs/nu2)**alpha[j]
        ax.plot(freqs, F1, linestyle='-', label=r'$\alpha = $'+str(alpha[j]))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel('aperture flux (Jy)')
    ax.set_xlabel('frequency (GHz)')
    #ax.legend()

    
def create_multi_plot(tab_path):
    tab = Table.read(tab_path)
    ind = get_ind(tab)

    cols = 6
    rows = len(ind)/6
    if rows > int(rows):
        rows = int(rows) + 1

    fig, axs = plt.subplots(rows, cols, sharex='all', figsize=(10,10))

    for i in range(len(ind)):
        ax = axs.flatten()[i]
        srcind = ind[i]
        plot_SED(ax, tab[srcind])

    plt.savefig('/home/jotter/nrao/plots/multi_SED_plot.png',dpi=300)
    
create_multi_plot('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_jun20.fits')

