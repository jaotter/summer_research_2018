#script to look at sources missed by Eisner 2018
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')

#Eisner FOV limits
eis_RA_upper = 83.836
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
srcI_ind = np.where(data['D_ID'][eis_ind] == 10)[0]
BN_ind = np.where(data['D_ID'][eis_ind] == 20)[0]
eis_ind = np.delete(eis_ind, np.array([srcI_ind, BN_ind]))

#load in eisner data:
eisner_data = Table.read('/users/jotter/summer_research_2018/tables/eisner_tbl.txt', format='ascii')
eis_coord_tab = Table.read('/users/jotter/summer_research_2018/tables/eis_coord_table.fits')

B3_coord = SkyCoord(ra=data['RA_B3']*u.deg, dec=data['DEC_B3']*u.deg)
eis_coord = SkyCoord(ra=eis_coord_tab['RA']*u.deg, dec=eis_coord_tab['DEC']*u.deg)

idx, d2d, d3d = eis_coord.match_to_catalog_sky(B3_coord)
match = np.where(d2d < 1*u.arcsec)[0] #match within 0.1 arcsec

data['Eis_ID'] = np.repeat(np.array('-',dtype='S8'), len(data))

for e_ind in match:
    data['Eis_ID'][idx[e_ind]] = eis_coord_tab['ID'][e_ind]

eis_match_ind = np.where(data['Eis_ID'] != '-')[0]

print(data['D_ID'][idx[match]])

fov_src_hist, bins = np.histogram(data['fwhm_maj_B3'][eis_ind], density=False)
matched_hist, b = np.histogram(data['fwhm_maj_B3'][eis_match_ind], bins, density=False)

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
ax1.set_xlabel('Deconvolved FWHM major (arcseconds)')

fov_src_hist_flux, bins = np.histogram(data['ap_flux_B3'][eis_ind], density=False)
matched_hist_flux, b = np.histogram(data['ap_flux_B3'][eis_match_ind], bins, density=False)

#fov_src_hist_flux, bins = np.histogram(np.log10(data['ap_flux_B3'][eis_ind]), density=False)
#matched_hist_flux, b = np.histogram(np.log10(data['ap_flux_B3'][eis_match_ind]), bins, density=False)


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

plt.savefig('plots/comb_hist_missed_srcs.png', dpi=400)



fig = plt.figure()

plt.errorbar(data['ap_flux_B3'][eis_ind], data['fwhm_maj_B3'][eis_ind], xerr=data['ap_flux_err_B3'][eis_ind], yerr=data['fwhm_maj_err_B3'][eis_ind], label='B3 sources in E18 FOV', linestyle='', marker='.')
plt.errorbar(data['ap_flux_B3'][eis_match_ind], data['fwhm_maj_B3'][eis_match_ind], xerr=data['ap_flux_err_B3'][eis_match_ind], yerr=data['fwhm_maj_err_B3'][eis_match_ind], label='B3 sources detected by E18', linestyle='', marker='.')
plt.legend()
plt.ylabel('Band 3 flux (Jy)')
plt.xlabel('Deconovlved FWHM major (as)')
plt.savefig('plots/size_flux_missed_srcs.png', dpi=400)
