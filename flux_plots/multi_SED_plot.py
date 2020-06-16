from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

    print(f'creating SED for source {srctab["D_ID"]}')
    
    fluxes = [srctab['ap_flux_B3']*1000, srctab['ap_flux_B6']*1000, srctab['ap_flux_B7']*1000]
    flux_errs = [srctab['ap_flux_err_B3'], srctab['ap_flux_err_B6'], srctab['ap_flux_err_B7']]

    F2 = fluxes[1]
    nu2 = freqs[1]
   
    ax.errorbar(freqs, fluxes, yerr=flux_errs, linestyle='', marker='o')
               
    for j in range(len(alpha)):
        F1 = F2*(freqs/nu2)**alpha[j]
        ax.plot(freqs, F1, linestyle='-', label=r'$\alpha = $'+str(alpha[j]))

    ax.text(0.1,0.82,srctab['D_ID'],transform=ax.transAxes,weight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.legend()

    
def create_multi_plot(tab_path, ncols=5):
    tab = Table.read(tab_path)
    ind = get_ind(tab)

    cols = ncols
    rows = len(ind)/ncols
    if rows > int(rows):
        rows = int(rows) + 1

    fig, axs = plt.subplots(rows, cols, sharex='all', figsize=(10,10))

    for i in range(len(axs.flatten())):
        ax = axs.flatten()[i]

        if i >= len(ind):
            ax.tick_params(axis='both',which='both',length=0,labelsize=0)
            
        else:
            srcind = ind[i]
            plot_SED(ax, tab[srcind])

            if len(axs.flatten()) - i <= ncols:
                ax.set_xlabel('Frequency (GHz)')
                ax.set_xticks([100,200,300])
                ax.set_xticklabels(['100', '200', '300'])
            if i % ncols == 0:
                ax.set_ylabel('Flux (mJy)')
                
    plt.subplots_adjust(wspace=0.35, hspace=0.3)
    plt.savefig('/home/jotter/nrao/plots/multi_SED_plot.pdf')

    
create_multi_plot('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_jun20.fits')

