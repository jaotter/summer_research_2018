from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.stats import ks_2samp
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')

irtab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
ir_did = irtab['ID']
ir_ind = []
for did in ir_did:
    ir_ind.append(np.where(data['ID'] == did)[0][0])

nonir_ind = np.setdiff1d(np.arange(0,len(data)), ir_ind)
    
#new_srcs = [9,10,22,24,36,37,39,45,55,59,61,71,72,79,81,84]
new_srcs = [8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124]

new_ind=[]
for src in new_srcs:
    new_ind.append(np.where(data['ID'] == src)[0][0])
new_ind=np.array(new_ind)

newtab = data[new_ind]
newtab.write('/home/jotter/nrao/summer_research_2018/tables/r0.5_new_det_may21.fits', overwrite=True)


#srcI_ind = np.where(data['D_ID'][eis_ind] == 30)[0]
#BN_ind = np.where(data['D_ID'][eis_ind] == 43)[0]

all_src_hist, bins = np.histogram(data['fwhm_maj_B3'], density=False)
new_src_hist, b = np.histogram(data['fwhm_maj_B3'][new_ind], bins=bins, density=False)


print('KS test for size:')
Dval, pval = ks_2samp(data['fwhm_maj_B3'][new_ind], data['fwhm_maj_B3'])
print('2 sample KS statistic: %f, p value: %f' % (Dval, pval))


#plotpts = []
#widths = []
#for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
#    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
#    widths.append((bins[b+1]-bins[b]))
                            
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(8,6))
fig.subplots_adjust(hspace=0.3)

widths = bins[1] - bins[0]

ax1.bar(bins[:-1], new_src_hist, widths, label='New sources', alpha=0.5, edgecolor='k', align='edge')
ax1.bar(bins[:-1], all_src_hist, widths, label='All sources', alpha=0.5, edgecolor='k', align='edge')
ax1.legend()
ax1.set_ylabel('Number')
ax1.set_xlabel('Non-deconvolved FWHM major (arcseconds)')

all_hist_flux, bins = np.histogram(np.log10(data['ap_flux_B3']), density=False)
new_src_hist_flux, b = np.histogram(np.log10(data['ap_flux_B3'][new_ind]), bins=bins, density=False)

print('KS test for flux:')
Dval, pval = ks_2samp(np.log10(data['ap_flux_B3'][new_ind]), np.log10(data['ap_flux_B3']))
print('2 sample KS statistic: %f, p value: %f' % (Dval, pval))
    
#plotpts = []
#widths = []
#for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
#    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
#    widths.append((bins[b+1]-bins[b]))

widths = bins[1] - bins[0]

ax2.bar(bins[:-1], new_src_hist_flux, widths, label='New sources', alpha=0.5, edgecolor='k', align='edge')
ax2.bar(bins[:-1], all_hist_flux, widths, label='All sources', alpha=0.5, edgecolor='k', align='edge')
ax2.legend()
ax2.set_ylabel('Number')
ax2.set_xlabel('log(Band 3 flux (Jy))')

#plt.savefig('/home/jotter/nrao/plots/comb_hist_new_srcs.pdf', dpi=300, bbox_inches='tight')



fig = plt.figure()

plt.errorbar(data['ap_flux_B3'], data['fwhm_maj_B3'], xerr=data['ap_flux_err_B3'], yerr=data['fwhm_maj_err_B3'], label='All sources', linestyle='', marker='.')
plt.errorbar(data['ap_flux_B3'][new_ind], data['fwhm_maj_B3'][new_ind], xerr=data['ap_flux_err_B3'][new_ind], yerr=data['fwhm_maj_err_B3'][new_ind], label='New sources', linestyle='', marker='.')
plt.legend(loc='upper left')
plt.loglog()
plt.xlabel('Band 3 flux (Jy)')
plt.ylabel('Non-deconovlved FWHM major (arcsecond)')
plt.savefig('/home/jotter/nrao/plots/size_flux_new_srcs.pdf', dpi=300, bbox_inches='tight')

fig = plt.figure()

plt.errorbar(data['ap_flux_B3'][nonir_ind], data['fwhm_maj_B3'][nonir_ind], xerr=data['ap_flux_err_B3'][nonir_ind], yerr=data['fwhm_maj_err_B3'][nonir_ind], label='OMC1', linestyle='', marker='.')
plt.errorbar(data['ap_flux_B3'][ir_ind], data['fwhm_maj_B3'][ir_ind], xerr=data['ap_flux_err_B3'][ir_ind], yerr=data['fwhm_maj_err_B3'][ir_ind], label='ONC', linestyle='', marker='.')
plt.legend(loc='upper left')
plt.loglog()
plt.xlabel('Band 3 flux (Jy)')
plt.ylabel('Non-deconovlved FWHM major (as)')
#plt.savefig('/home/jotter/nrao/plots/size_flux_onc_omc1.pdf', dpi=300, bbox_inches='tight')
