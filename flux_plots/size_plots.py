from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u

import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

from lifelines import KaplanMeierFitter

def dust_mass_KM():
    vt19_data = Table.read('/home/jotter/nrao/tables/VT19.txt', format='latex')
    eis_data = Table.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='ascii')
    dmass_data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_apr20_calc_vals.fits')

    vt19_dmass_raw = vt19_data['Mass']
    eis_dmass_raw = eis_data['M_dust^a']
    B3_dmass = np.log10(dmass_data['dust_mass_B3'])
    B7_dmass = np.log10(dmass_data['dust_mass_B7'])
    B7_dmass = B7_dmass[np.where(np.isnan(B7_dmass) == False)[0]]

    vt19_dmass = []
    for dm in vt19_dmass_raw:
        vt19_dmass.append(dm.split()[0][1:])

    vt19_dmass = np.log10(np.array(vt19_dmass[1:],dtype='float'))
    
    #eis_dmass = []
    #for dm in eis_dmass_raw:
    #    eis_dmass.append(dm.split()[0])
    #eis_dmass = np.log10(np.array(eis_dmass, dtype='float'))
    #eis_dmass = eis_dmass[np.where(np.isinf(eis_dmass) == False)[0]]

    kmf = KaplanMeierFitter()
    kmf.fit(vt19_dmass)

    fig = plt.figure()
    ax = kmf.plot()
    plt.savefig('/home/jotter/nrao/plots/VT19_KM_plot.png')
    
def KM_hist(hist):
    #hist is histogram to compute KM estimator for (deaths), bins is right edge of hist bins (time)
    Si = []
    Sprob = []
    surv = np.sum(hist)
    for i in range(len(hist)):
        val = 1 - hist[i]/surv
        surv = surv - hist[i]
        print(surv, val, hist[i])
        Si.append(val)
        Sprob.append(np.prod(Si))

    return Sprob
    
    
def dust_mass_comp():
    vt19_data = Table.read('/home/jotter/nrao/tables/VT19.txt', format='latex')
    eis_data = Table.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='ascii')
    dmass_data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_apr20_calc_vals.fits')

    vt19_dmass_raw = vt19_data['Mass']
    eis_dmass_raw = eis_data['M_dust^a']
    B3_dmass = np.log10(dmass_data['dust_mass_B3'])
    B7_dmass = np.log10(dmass_data['dust_mass_B7'])
    B7_dmass = B7_dmass[np.where(np.isnan(B7_dmass) == False)[0]]

    vt19_dmass = []
    for dm in vt19_dmass_raw:
        vt19_dmass.append(dm.split()[0][1:])

    vt19_dmass = np.log10(np.array(vt19_dmass[1:],dtype='float'))
    
    eis_dmass = []
    for dm in eis_dmass_raw:
        eis_dmass.append(dm.split()[0])
    eis_dmass = np.log10(np.array(eis_dmass, dtype='float'))
    eis_dmass = eis_dmass[np.where(np.isinf(eis_dmass) == False)[0]]
    
    print(vt19_dmass, eis_dmass, B3_dmass)
    all_masses_B3 = np.concatenate((vt19_dmass, eis_dmass, B3_dmass))
    all_hist_B3, bins_B3 = np.histogram(all_masses_B3)

    vt_hist_B3, b = np.histogram(vt19_dmass, bins=bins_B3)
    eis_hist_B3, b = np.histogram(eis_dmass, bins=bins_B3)
    B3_hist, b = np.histogram(B3_dmass, bins=bins_B3)
    eis_B3_combined_hist = eis_hist_B3+B3_hist
    
    plotpts_B3 = []
    widths_B3 = []
    for b in range(len(bins_B3[:-1])): #creating points to plot - midpoints of bins
        plotpts_B3.append(bins_B3[b] + (bins_B3[b+1]-bins_B3[b])/2)
        widths_B3.append((bins_B3[b+1]-bins_B3[b]))

    all_masses_B7 = np.concatenate((vt19_dmass, eis_dmass, B7_dmass))
    all_hist_B7, bins_B7 = np.histogram(all_masses_B7)

    vt_hist_B7, b = np.histogram(vt19_dmass, bins=bins_B7)
    eis_hist_B7, b = np.histogram(eis_dmass, bins=bins_B7)
    B7_hist, b = np.histogram(B7_dmass, bins=bins_B7)
    eis_B7_combined_hist = eis_hist_B7+B7_hist
    
    plotpts_B7 = []
    widths_B7 = []
    for b in range(len(bins_B7[:-1])): #creating points to plot - midpoints of bins
        plotpts_B7.append(bins_B7[b] + (bins_B7[b+1]-bins_B7[b])/2)
        widths_B7.append((bins_B7[b+1]-bins_B7[b]))

    vtdist = KM_hist(vt_hist_B3)
    eisdist = KM_hist(eis_hist_B3)
    B3dist = KM_hist(B3_hist)
    eisB3dist = KM_hist(eis_B3_combined_hist)

    fig = plt.figure()
    plt.plot(bins_B3[1:], vtdist, linestyle='--', label='VT19', marker='o')
    plt.plot(bins_B3[1:], eisdist, linestyle=':', label='E18', marker='o')
    plt.plot(bins_B3[1:], B3dist, linestyle='-', label='band 3', marker='o')
    plt.xlabel(r'$\log(M_{dust}/M_\oplus)$')
    plt.ylabel('Probability of greater mass')
    plt.legend()
    plt.savefig('/home/jotter/nrao/plots/KM_dust_mass_B3.png', dpi=400)

    fig = plt.figure()
    plt.plot(bins_B3[1:], vtdist, linestyle='--', label='VT19', marker='o')
    plt.plot(bins_B3[1:], eisB3dist, linestyle='-', label='band 3 and E18', marker='o')
    plt.xlabel(r'$\log(M_{dust}/M_\oplus)$')
    plt.ylabel('Probability of greater mass')
    plt.legend()
    plt.savefig('/home/jotter/nrao/plots/KM_dust_mass_B3_combined.png', dpi=400)

    
    fig = plt.figure()
    plt.bar(plotpts_B3, vt_hist_B3, widths_B3, edgecolor = 'black', label='VT19', alpha=0.5)
    plt.bar(plotpts_B3, eis_hist_B3, widths_B3, edgecolor = 'black', label='E18', alpha=0.5)
    plt.bar(plotpts_B3, B3_hist, widths_B3, edgecolor = 'black', label='band 3', alpha=0.5)
    plt.xlabel(r'$\log(M_{dust}/M_\oplus)$')
    plt.ylabel('number of disks')
    plt.legend()
    plt.savefig('/home/jotter/nrao/plots/dust_mass_hist_B3.png', dpi=400)

    '''fig = plt.figure()
    plt.bar(plotpts_B7, vt_hist_B7, widths_B7, edgecolor = 'black', label='VT19', alpha=0.5)
    plt.bar(plotpts_B7, eis_hist_B7, widths_B7, edgecolor = 'black', label='E18', alpha=0.5)
    plt.bar(plotpts_B7, B7_hist, widths_B7, edgecolor = 'black', label='band 7', alpha=0.5)
    plt.xlabel(r'$\log(M_{dust}/M_\oplus)$')
    plt.ylabel('number of disks')
    plt.legend()
    plt.savefig('/home/jotter/nrao/plots/dust_mass_hist_B7.png', dpi=400)'''
    
    fig = plt.figure()
    plt.bar(plotpts_B3, vt_hist_B3, widths_B3, edgecolor = 'black', label='VT19', alpha=0.5)
    plt.bar(plotpts_B3, eis_B3_combined_hist, widths_B3, edgecolor = 'black', label='E18 and band 3', alpha=0.5)
    plt.xlabel(r'$\log(M_{dust}/M_\oplus)$')
    plt.ylabel('number of disks')
    plt.legend()
    plt.savefig('/home/jotter/nrao/plots/dust_mass_hist_B3_combined.png', dpi=400)

    
def disk_size_hist(arrs, labels, filename):
    fig = plt.figure()
    for a in range(len(arrs)):
        size_arr = arrs[a]
        size_arr = size_arr[np.isnan(size_arr)==False]
        print(len(size_arr))
        hist, bins = np.histogram(size_arr, density=False)
        plotpts = []
        widths = []
        for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
            plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
            widths.append((bins[b+1]-bins[b]))
    
        plt.bar(plotpts, hist, widths, edgecolor = 'black', label=labels[a], alpha=0.5)
    plt.legend()
    plt.xlabel('deconvolved FWHM major (as)')
    plt.ylabel('number of disks')
    plt.xlim(0,0.35)
    plt.savefig('/home/jotter/nrao/plots/size_plots/'+filename, dpi=500)

def disk_size_hist_3panel(arrs, xlabels, filename, nbins=10):
    f, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(12,12), sharex='col')
    ax = [ax1, ax2, ax3]

    d = 414*u.pc
    d = d.to(u.AU)
    
    total = np.concatenate(arrs)
    total = total[np.logical_not(np.isnan(total))]
    total_hist, bins = np.histogram(total, bins=nbins)

    B3ind = np.where(np.isnan(arrs[0]) == False)
    B6ind = np.where(np.isnan(arrs[1]) == False)
    B7ind = np.where(np.isnan(arrs[2]) == False)
    allind = reduce(np.intersect1d, (B3ind, B6ind, B7ind))
    print(len(allind))
    plotpts = []
    widths = []
    for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
        plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
        widths.append((bins[b+1]-bins[b]))
    
    for a in range(len(arrs)):
        size_arr = arrs[a]
        size_arr = size_arr[np.isnan(size_arr)==False]
        hist, b = np.histogram(size_arr, bins, density=False)
        allhist, b = np.histogram(arrs[a][allind], bins, density=False)

        ax[a].bar(plotpts, hist, widths, edgecolor = 'black', alpha=0.5, label=f'{xlabels[a]} sources')
        ax[a].bar(plotpts, allhist, widths, edgecolor='black', alpha=0.5, label=f'{xlabels[a]} sources detected in all bands')
        ax[a].set_xlim(0,0.27)
        ax[a].set_ylim(0,15)
        ax[a].set_ylabel('number of disks', fontsize=18)
        ax[a].legend(fontsize=16)
        ax[a].yaxis.set_major_locator(plt.MaxNLocator(5))
        
    altax = ax1.twiny()
    xlim = 0.27*u.arcsec.to(u.rad)
    xlim *= d
    altax.set_xlim(0, xlim.value)
    altax.set_xlabel('deconvolved FWHM major (AU)', fontsize=18)
    ax3.set_xlabel('deconvolved FWHM major (as)', fontsize=18)

    ax1.set_yticks(ax1.get_yticks()[1:])
    ax2.set_yticks(ax2.get_yticks()[1:])

    altax.tick_params(axis='x', which='both', labelsize='large')
    ax1.tick_params(axis='y', which='both', labelsize='large')
    ax2.tick_params(axis='y', which='both', labelsize='large')
    ax3.tick_params(axis='both', which='both', labelsize='large')
    
    plt.tight_layout()
    f.subplots_adjust(hspace=0)

    plt.savefig('/home/jotter/nrao/plots/size_plots/'+filename, dpi=500)


def R_hist_eisner(size_arr, xlabel, label, filename, nbins=10):
    fig = plt.figure()

    d = 414*u.pc
    d = d.to(u.AU)
    
    size_arr = size_arr[np.isnan(size_arr)==False]*u.arcsec
    size_arr = (size_arr.to(u.rad)*d).value
    hist, bins = np.histogram(size_arr, density=False, bins=nbins)
    plotpts = []
    widths = []
    for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
        plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
        widths.append((bins[b+1]-bins[b]))

    eisner_data = ascii.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='tab')
    eisner_ind = np.where(eisner_data['R_disk'] != '<5')[0]
    eisner_R = [float(x.split()[0])*2 for x in eisner_data['R_disk'][eisner_ind]]
    eis_hist, bins = np.histogram(eisner_R, density=False, bins=bins)
    
    plt.bar(plotpts, hist, widths, edgecolor = 'black', label=label, alpha=0.4)
    plt.bar(plotpts, eis_hist, widths, edgecolor='black', label='Eisner sizes', alpha=0.4)
    plt.legend()
    plt.xlabel('deconvolved FWHM major (AU)')
    plt.ylabel('number of disks')
    #plt.xlim(0,0.35)
    #plt.style.use(mpl_style.style1)
    plt.savefig('/home/jotter/nrao/plots/size_plots/'+filename, dpi=400)
    
def size_comp(conv_arrs, deconv_arrs, conv_errs, deconv_errs, labels, filename):
    #plot sizes of disks in two bands, arr1 and arr2 should have same length and sources
    #labels - arr1 label, then arr2

    fig = plt.figure()

    deconv_inds = np.where(np.isnan(deconv_arrs[0]) == False)[0]

    for ind in deconv_inds:
        plt.errorbar([deconv_arrs[0][ind], conv_arrs[0][ind]], [deconv_arrs[1][ind], conv_arrs[1][ind]], xerr=[deconv_errs[0][ind], conv_errs[0][ind]], yerr=[deconv_errs[1][ind], conv_errs[1][ind]], marker='o', color='r')

    plt.errorbar(conv_arrs[0], conv_arrs[1], xerr=conv_errs[0], yerr=conv_errs[1], marker='o',linestyle='', label='convolved sizes')
    
        
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.xlim(0.0,0.3)
    plt.ylim(0.0,0.3)
    plt.plot(np.arange(0,1,0.1), np.arange(0,1,0.1), color='k')
    plt.legend()
    #plt.style.use(mpl_style.style1)
    plt.savefig('/home/jotter/nrao/plots/size_plots/'+filename, dpi=400)

def size_comp_simple(arrs, errs, labels, filename):
    
    #plot sizes of disks in two bands, arr1 and arr2 should have same length and sources
    #labels - arr1 label, then arr2

    fig = plt.figure()

    ind1 = np.where(arrs[0]-3*errs[0] > arrs[1]+3*errs[1])[0]
    ind2 = np.where(arrs[0]+3*errs[0] < arrs[1]-3*errs[1])[0]
    ind = np.concatenate((ind1, ind2))
    print(ind)
    plt.errorbar(arrs[0], arrs[1], xerr=3*errs[0], yerr=3*errs[1], linestyle='', marker='.')    
    plt.errorbar(arrs[0][ind], arrs[1][ind], xerr=3*errs[0][ind], yerr=3*errs[1][ind], linestyle='', marker='.')    
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.xlim(0.0,0.3)
    plt.ylim(0.0,0.3)
    plt.plot(np.arange(0,1,0.1), np.arange(0,1,0.1), color='k')
    #plt.style.use(mpl_style.style1)
    plt.savefig('/home/jotter/nrao/plots/size_plots/'+filename, dpi=400)

def size_comp_eisner(filename):
    fig = plt.figure()

    data = Table.read('/home/jotter/nrao/tables/eis_r0.5_apr20_match.fits')
    
    eisner_data = ascii.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='tab')
    eis_UL_ind = np.where(eisner_data['R_disk'] == '<5')[0]
    eisner_R_str = eisner_data['R_disk']
    eisner_R_str[eis_UL_ind] = '0.0'
    eisner_R = [float(x.split()[0])*2 for x in eisner_R_str]
    eisner_R_err = [float(x.split()[-1])*2 for x in eisner_R_str]
    
    data_R = data['fwhm_maj_deconv_B3']*u.arcsec
    data_R = (data_R.to(u.rad)*(414*u.pc).to(u.au)).value
    data_R_err = data['fwhm_maj_err_B3']*u.arcsec
    data_R_err = (data_R_err.to(u.rad)*(414*u.pc).to(u.au)).value

    deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0]
    data_R = data_R[deconv_ind]
    data_R_err = data_R_err[deconv_ind]
    
    match_eis_ind = []
    for ID in data['ID'][deconv_ind]:
        match_eis_ind.append(np.where(eisner_data['ID'] == ID.strip())[0][0])
        print(ID, match_eis_ind[-1])
        
    eisner_R = np.array(eisner_R)[match_eis_ind]
    eisner_R_err = np.array(eisner_R_err)[match_eis_ind]

    ind1 = np.where(data_R+3*data_R_err < eisner_R - 3*eisner_R_err)[0]
    ind2 = np.where(data_R-3*data_R_err > eisner_R + 3*eisner_R_err)[0]
    ind = np.concatenate((ind1, ind2))
    
    plt.errorbar(data_R, eisner_R, xerr=3*data_R_err, yerr=3*eisner_R_err, linestyle='', marker='.')
    plt.errorbar(data_R[ind], eisner_R[ind], xerr=3*data_R_err[ind], yerr=3*eisner_R_err[ind], linestyle='', marker='.')
    plt.xlabel('B3 FWHM (AU)')
    plt.ylabel('Eisner et al 2018 FWHM (AU)')
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.plot(np.arange(0,101,50), np.arange(0,101,50), color='k')
    #plt.style.use(mpl_style.style1)
    plt.savefig('/home/jotter/nrao/plots/size_plots/'+filename, dpi=400)
    
def img_size_comp(arrs, err_arrs,  labels):
    fig = plt.figure()
    inds1 = np.where(err_arrs[0] < 0.2)[0]
    inds2 = np.where(err_arrs[1] < 0.2)[0]
    inds = np.intersect1d(inds1, inds2)
    plt.errorbar(arrs[0][inds], arrs[1][inds], xerr=err_arrs[0][inds], yerr=err_arrs[1][inds], linestyle='', marker='o')
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
