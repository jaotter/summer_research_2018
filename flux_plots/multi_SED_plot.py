from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.table import Table
from scipy.optimize import curve_fit
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from functools import reduce
from astropy.coordinates import Angle, SkyCoord

from matplotlib.ticker import ScalarFormatter

def f(x, m, b):
    return m*x + b

def get_ind(table): #returns indices where sources are detected in all bands
    b6_ind = np.where(np.isnan(table['ap_flux_B6']) == False)[0]
    b6_ulim = np.where(np.isnan(table['B6_flux_ulim']) == False)[0]
    b7_ind = np.where(np.isnan(table['ap_flux_B7']) == False)[0]
    b7_ulim = np.where(np.isnan(table['B7_flux_ulim']) == False)[0]

    b6 = np.concatenate((b6_ind, b6_ulim))
    b7 = np.concatenate((b7_ind, b7_ulim))

    both_ulim = np.intersect1d(b6_ulim, b7_ulim)
    
    ind = np.intersect1d(b6, b7)
    ind = np.concatenate((ind, [4, 6, 11, 3, 19]))
    ind = np.sort(ind)

    return ind


def plot_SED(ax, srcid, fluxes, flux_errs, ulims, e18_b7=False):
    freqs = np.array([98, 223.5, 339.7672758867]) #*u.GHz

    print(f'creating SED for source {srcid}')
    
    flux_errs = np.log10(np.exp(1)) * flux_errs / fluxes #error propagating for log
    fluxes = np.log10(fluxes)
    ulims = np.log10(ulims)
    fluxind = np.where(np.isnan(fluxes) == False)[0]
    ulimind = np.where(np.isnan(ulims) == False)[0]

    ax.errorbar(freqs[fluxind], fluxes[fluxind], yerr=flux_errs[fluxind], linestyle='', marker='o', markersize=2, color='tab:blue')
    
    if e18_b7 == True:
        ax.errorbar(freqs[-1], fluxes[-1], yerr=flux_errs[-1], linestyle='', marker='s', markersize=3, color='tab:green')
        
    if len(ulimind) > 0:
        ax.plot(freqs[ulimind], ulims[ulimind], linestyle='', marker='v', color='tab:orange', markersize=2)

    if len(fluxind) > 1:
        popt, pcov = curve_fit(f, np.log10(freqs[fluxind]), fluxes[fluxind], sigma=flux_errs[fluxind])
        ax.plot(freqs, np.log10(freqs)*popt[0] + popt[1], linestyle='-', color='k', alpha=0.5)
        perr = np.sqrt(np.diag(pcov))

        if srcid == 14 or srcid == 3 or srcid == 19:
            freqdiff = np.log10(freqs[2]) - np.log10(freqs[0])
            perr[0] = np.sqrt((flux_errs[0]/(freqdiff*np.log(10)*fluxes[0]))**2 + (flux_errs[2]/(freqdiff*np.log(10)*fluxes[2]))**2)
            specind_str = str(np.round(popt[0],1))+r'$\pm$'+str(np.round(perr[0],1))
            if srcid == 14:
                ax.text(0.45,0.4,specind_str,transform=ax.transAxes, fontsize=9)
            else:
                ax.text(0.45,0.07,specind_str,transform=ax.transAxes, fontsize=9)
        else:
            specind_str = str(np.round(popt[0],1))+r'$\pm$'+str(np.round(perr[0],1))
            ax.text(0.45,0.07,specind_str,transform=ax.transAxes, fontsize=9)

    else:
        popt = (np.nan, np.nan)
        perr = (np.nan, np.nan)

    alpha=2
    F2 = fluxes[1]
    nu2 = freqs[1]
    if np.isnan(F2) == True:
        F2 = fluxes[0]
        nu2 = freqs[0]
    F1 = np.log10(np.power(10,F2)*(freqs/nu2)**alpha)
    ax.plot(freqs, F1, linestyle='--', color='tab:red', alpha=0.5)    
    ax.text(0.1,0.82,srcid,transform=ax.transAxes,weight='bold')

    ax.set_xscale('log')
    #ax.set_yscale('log')

    return popt, perr

    
def create_multi_plot(tab_path, ncols=6):
    tab = Table.read(tab_path)
    ind = get_ind(tab)

    fit_alpha = []
    fit_alpha_err = []
    
    cols = ncols
    rows = len(ind)/ncols
    if rows > int(rows):
        rows = int(rows) + 1
    
    fig, axs = plt.subplots(int(rows), cols, figsize=(12,10), sharex='all')

    #manually adding E18 fluxes:
    e18_fluxes = np.array([11.2, 6.7, 13, 6.6, 1.6])
    e18_flux_errs = np.array([0.4, 0.7, 1.3, 0.6, 0.3])
    e18_srcids = np.array([3, 4, 6, 11, 19])

    for i in range(len(axs.flatten())):
        ax = axs.flatten()[i]

        if i >= len(ind):
            #ax.tick_params(axis='both',which='both',length=0,labelsize=0)
            ax.set_axis_off()
        else:   
            srcind = ind[i]
            fluxes = np.array([tab['ap_flux_B3'][srcind], tab['ap_flux_B6'][srcind], tab['ap_flux_B7'][srcind]]) * 1000
            flux_errs = np.array([tab['ap_flux_err_B3'][srcind], tab['ap_flux_err_B6'][srcind], tab['ap_flux_err_B7'][srcind]]) * 1000
            sys_flux_errs = np.sqrt(np.array(flux_errs)**2 + (0.2*np.array(fluxes))**2)
            ulims = np.array([np.nan, tab['B6_flux_ulim'][srcind], tab['B7_flux_ulim'][srcind]])

            if srcind in e18_srcids:
                e18ind = np.where(e18_srcids == srcind)[0]
                fluxes[-1] = e18_fluxes[e18ind]
                sys_flux_errs[-1] = e18_flux_errs[e18ind]
                fit_vals, fit_errs = plot_SED(ax, tab['ID'][srcind], fluxes, sys_flux_errs, ulims, e18_b7=True)
            else:
                fit_vals, fit_errs = plot_SED(ax, tab['ID'][srcind], fluxes, sys_flux_errs, ulims, e18_b7=False)

            fit_alpha.append(fit_vals[0])
            fit_alpha_err.append(fit_errs[0])

            #ax.ticklabel_format(style='plain')
            #ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
            #ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
            #ax.yaxis.set_minor_formatter(NullFormatter())
            #ax.yaxis.set_minor_formatter(StrMethodFormatter('{x:.0f}'))

            #ax.xaxis.set_major_formatter(ScalarFormatter())
            #ax.yaxis.set_major_formatter(ScalarFormatter())
            
            if len(axs.flatten()) - i <= ncols:
                ax.set_xticks([100,200,300])
                ax.set_xticklabels(['100', '200', '300'])
                
    fig.text(0.5, 0.06, 'Frequency (GHz)', ha='center', fontsize=14)
    fig.text(0.06, 0.5, 'log Flux (mJy)', va='center', rotation='vertical', fontsize=14)
    
    plt.subplots_adjust(wspace=0.4, hspace=0.3)
    plt.savefig('/home/jotter/nrao/plots/multi_SED_plot_E18.pdf', bbox_inches='tight')


    calc_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_may21_calc_vals_mask.fits')
    calc_tab['alpha_fit'] = np.repeat(np.nan, len(calc_tab))
    calc_tab['alpha_fit'][ind] = fit_alpha
    calc_tab['alpha_fit_err'] = np.repeat(np.nan, len(calc_tab))
    calc_tab['alpha_fit_err'][ind] = fit_alpha_err    

    calc_tab.write('/home/jotter/nrao/summer_research_2018/tables/r0.5_may21_calc_vals_mask_alpha.fits', overwrite=True)
    
create_multi_plot('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')

