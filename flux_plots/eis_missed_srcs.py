#script to look at sources missed by Eisner 2018
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.stats import ks_2samp
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_jun20.fits')

#Eisner FOV limits

'''eis_RA_upper = 83.836
eis_RA_lower = 83.805
eis_Dec_upper = -5.3715
eis_Dec_lower = -5.402

RA_ind1 = np.where(data['RA_B3'] > eis_RA_lower)[0]
RA_ind2 = np.where(data['RA_B3'] < eis_RA_upper)[0]
RA_ind = np.intersect1d(RA_ind1,RA_ind2)

Dec_ind1 = np.where(data['DEC_B3'] > eis_Dec_lower)[0]
Dec_ind2 = np.where(data['DEC_B3'] < eis_Dec_upper)[0]
Dec_ind = np.intersect1d(Dec_ind1,Dec_ind2)

eis_ind = np.intersect1d(RA_ind, Dec_ind)
'''

eis_fov_srcs = [47,44,21,19,7,4,73,2,5,3,6,8,13,24,48,14,17]
eis_ind=[]
for src in eis_fov_srcs:
    eis_ind.append(np.where(data['D_ID'] == src)[0][0])
eis_ind=np.array(eis_ind)

srcI_ind = np.where(data['D_ID'][eis_ind] == 30)[0]
BN_ind = np.where(data['D_ID'][eis_ind] == 43)[0]
eis_ind = np.delete(eis_ind, np.array([srcI_ind, BN_ind]))
eis_ind = np.delete(eis_ind, [17])

#load in eisner data:
eisner_data = Table.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='ascii')
eis_coord_tab = Table.read('/home/jotter/nrao/tables/eis_coord_table.fits')

eis_match_tab = Table.read('/home/jotter/nrao/tables/eis_r0.5_jun20_match.fits')

fov_src_hist, bins = np.histogram(data['fwhm_maj_B3'][eis_ind], density=False)
matched_hist, b = np.histogram(eis_match_tab['fwhm_maj_B3'], bins, density=False)

print('KS test for size:')
Dval, pval = ks_2samp(data['fwhm_maj_B3'][eis_ind], eis_match_tab['fwhm_maj_B3'])
print('2 sample KS statistic: %f, p value: %f' % (Dval, pval))


plotpts = []
widths = []
for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))
                            
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(8,6))
fig.subplots_adjust(hspace=0.3)


ax1.bar(plotpts, fov_src_hist, widths, label='B3 sources in E18 FOV', alpha=0.5, edgecolor='k')
ax1.bar(plotpts, matched_hist, widths, label='B3 sources detected by E18', alpha=0.5, edgecolor='k')
ax1.legend()
ax1.set_ylabel('Number')
ax1.set_xlabel('Non-deconvolved FWHM major (arcseconds)')

fov_src_hist_flux, bins = np.histogram(data['ap_flux_B3'][eis_ind], density=False)
matched_hist_flux, b = np.histogram(eis_match_tab['ap_flux_B3'], bins, density=False)

#fov_src_hist_flux, bins = np.histogram(np.log10(data['ap_flux_B3'][eis_ind]), density=False)
#matched_hist_flux, b = np.histogram(np.log10(data['ap_flux_B3'][eis_match_ind]), bins, density=False)

print('KS test for flux:')
Dval, pval = ks_2samp(data['ap_flux_B3'][eis_ind], eis_match_tab['ap_flux_B3'])
print('2 sample KS statistic: %f, p value: %f' % (Dval, pval))

maxflux_ind = np.where(data['ap_flux_B3'] == np.max(data['ap_flux_B3'][eis_ind]))
print(f'brightest source: {data["D_ID"][maxflux_ind]}')

    
plotpts = []
widths = []
for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))
                            
ax2.bar(plotpts, fov_src_hist_flux, widths, label='B3 sources in E18 FOV', alpha=0.5, edgecolor='k')
ax2.bar(plotpts, matched_hist_flux, widths, label='B3 sources detected by E18', alpha=0.5, edgecolor='k')
ax2.legend()
ax2.set_ylabel('Number')
ax2.set_xlabel('Band 3 flux (Jy)')

plt.savefig('/home/jotter/nrao/plots/E18_comp/comb_hist_missed_srcs2.pdf', dpi=400)



fig = plt.figure()

plt.errorbar(data['ap_flux_B3'][eis_ind], data['fwhm_maj_B3'][eis_ind], xerr=data['ap_flux_err_B3'][eis_ind], yerr=data['fwhm_maj_err_B3'][eis_ind], label='B3 sources in E18 FOV', linestyle='', marker='.')
plt.errorbar(eis_match_tab['ap_flux_B3'], eis_match_tab['fwhm_maj_B3'], xerr=eis_match_tab['ap_flux_err_B3'], yerr=eis_match_tab['fwhm_maj_err_B3'], label='B3 sources detected by E18', linestyle='', marker='.')
plt.legend()
plt.xlabel('Band 3 flux (Jy)')
plt.ylabel('Non-deconovlved FWHM major (as)')
plt.savefig('/home/jotter/nrao/plots/E18_comp/size_flux_missed_srcs2.pdf', dpi=400)
