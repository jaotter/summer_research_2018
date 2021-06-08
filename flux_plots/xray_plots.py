import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt


tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
tab_ir = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')

coup = Table.read('/home/jotter/nrao/summer_research_2018/tables/COUP_r0.5_may21.fits')
#coup_ir = Table.read('/home/jotter/nrao/summer_research_2018/tables/COUP_r0.5_may21_IRmatch.fits')

fb = Table.read('/home/jotter/nrao/summer_research_2018/tables/Forbrich2016_r0.5_may21.fits')
#fb_ir = Table.read('/home/jotter/nrao/summer_research_2018/tables/Forbrich2016_r0.5_jun20_IRmatch.fits')

'''
tab_ir_ind = []
for ir_did in tab_ir['D_ID']:
    tab_ir_ind.append(np.where(tab['D_ID'] == ir_did)[0][0])

tab_nonir_ind = np.delete(np.arange(0,len(tab)), tab_ir_ind)
tab_nonir = tab[tab_nonir_ind]


coup_nonir_dids = [18,29,32,33,35,41,45,54,58,81,83]
coup_nonir_ind = []
for coup_nonir in coup_nonir_dids:
    coup_nonir_ind.append(np.where(tab_nonir['D_ID'] == coup_nonir)[0][0])

fb_nonir_dids = [19, 28, 30, 31, 32, 41, 54, 56, 58, 67, 80]
fb_nonir_ind = []
for fb_nonir in fb_nonir_dids:
    fb_nonir_ind.append(np.where(tab_nonir['D_ID'] == fb_nonir)[0][0])
    
both_nonir_ind = np.array(set(np.concatenate((coup_nonir_ind, fb_nonir_ind))))
'''
coup_dids = coup['ID']
coup_ind = []
for coup_did in coup_dids:
    coup_ind.append(np.where(tab['ID'] == coup_did)[0][0])

fb_dids = fb['ID']
fb_ind = []
for fb_did in fb_dids:
    fb_ind.append(np.where(tab['ID'] == fb_did)[0][0])
    
''' 
##size hist

nonan_ind = np.where(np.isnan(tab_nonir['fwhm_maj_deconv_B3']) == False)[0]

print(tab_nonir['fwhm_maj_deconv_B3'])
print(tab_nonir['fwhm_maj_deconv_B3'][nonan_ind])
all_hist, bins = np.histogram(tab_nonir['fwhm_maj_deconv_B3'][nonan_ind])
print(bins)

plt.figure()
plt.hist(tab_nonir['fwhm_maj_deconv_B3'], bins=bins, label='OMC1', alpha=0.4)
plt.hist(tab_nonir['fwhm_maj_deconv_B3'][coup_nonir_ind], bins=bins, label='COUP', alpha=0.4)
plt.hist(tab_nonir['fwhm_maj_deconv_B3'][fb_nonir_ind], bins=bins, label='Forbrich2016', alpha=0.4)

plt.legend()
plt.savefig('/home/jotter/nrao/plots/OMC1_COUP_Fb_size.png')


##flux hist

#all_hist, bins = np.histogram(tab_nonir['ap_flux_B3'])

bins = np.linspace(0,0.01, 8)

plt.figure()
plt.hist(tab_nonir['ap_flux_B3'], bins=bins, label='OMC1', alpha=0.4)
plt.hist(tab_nonir['ap_flux_B3'][coup_nonir_ind], bins=bins, label='COUP', alpha=0.4)
plt.hist(tab_nonir['ap_flux_B3'][fb_nonir_ind], bins=bins, label='Forbrich2016', alpha=0.4)

plt.legend()
plt.savefig('/home/jotter/nrao/plots/OMC1_COUP_Fb_flux.png')

##all sources
'''
bins = np.linspace(0,0.01, 8)

plt.figure()
plt.hist(tab['ap_flux_B3'], bins=bins, label='all', alpha=0.4)
plt.hist(tab['ap_flux_B3'][coup_ind], bins=bins, label='COUP', alpha=0.4)
plt.hist(tab['ap_flux_B3'][fb_ind], bins=bins, label='Forbrich2016', alpha=0.4)

plt.legend()
plt.savefig('/home/jotter/nrao/plots/all_Coup_Fb_flux.png')
