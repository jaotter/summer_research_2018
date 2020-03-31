from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.table import Table
from astropy.nddata import Cutout2D
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

B3fl = fits.open('../../images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
header = B3fl[0].header
wcs = WCS(header).celestial
img = B3fl[0].data.squeeze()

IR_tab = Table.read('../../tables/A11_MLLA_r0.5_HC2000_LRY2000_matched_FOV.fits')

unmatched_IR_ind = np.where(np.isnan(IR_tab['RA_B3']) == True)[0]
matched_IR_ind = np.where(np.isnan(IR_tab['RA_B3']) == False)[0]

#first for unmatched IR ind -> get B3 upper lim
IR_um = IR_tab[unmatched_IR_ind]
um_pix_coords = wcs.all_world2pix(IR_um['RA']*u.degree, IR_um['DEC']*u.degree, 0)

img_inds = np.mgrid[0:len(img),0:len(img[0])]

rad_as = 0.4*u.arcsec
pix_scale = proj_plane_pixel_scales(wcs)
rad_pix = (rad_as.to(u.degree)/pix_scale[0]).value


B3_upper_lim = []
for ind in range(len(IR_um)):
    pix_x = int(um_pix_coords[0][ind])
    pix_y = int(um_pix_coords[1][ind])
    print(pix_x, pix_y)
    print(IR_um[ind]['MLLA'])

    img_cutout = Cutout2D(img, (pix_x,pix_y), 300)

    cutout_inds = np.mgrid[0:len(img_cutout.data),0:len(img_cutout.data[0])]
    circ_inds = np.where(np.sqrt((cutout_inds[0] - len(img_cutout.data)/2)**2 + (cutout_inds[1]-len(img_cutout.data[0])/2)**2) < rad_pix)


    std = np.nanstd(img_cutout.data[circ_inds])
    print(std)
    B3_upper_lim.append(3*std)
    
    
MLLA_k_ind = np.where(np.isnan(IR_tab['Kmag1'])==False)
HC2000_k_ind = np.where(np.isnan(IR_tab['Kmag_1'])==False)
LRY_k_ind = np.where(np.isnan(IR_tab['Kmag_2'])==False)

Kmag = np.repeat(np.nan, len(IR_tab))
Kmag_err = np.repeat(np.nan, len(IR_tab))

Kmag[LRY_k_ind] = IR_tab['Kmag_2'][LRY_k_ind]
Kmag_err[LRY_k_ind] = np.repeat(0.2, len(LRY_k_ind))
Kmag[HC2000_k_ind] = IR_tab['Kmag_1'][HC2000_k_ind]
Kmag_err[HC2000_k_ind] = IR_tab['e_Kmag'][HC2000_k_ind]
Kmag[MLLA_k_ind] = IR_tab['Kmag1'][MLLA_k_ind]
Kmag_err[MLLA_k_ind] = IR_tab['e_Kmag1'][MLLA_k_ind]

IR_m = IR_tab[matched_IR_ind]
Kmag_m = Kmag[matched_IR_ind]
Kmag_m_err = Kmag_err[matched_IR_ind]

Kmag_um = Kmag[unmatched_IR_ind]
Kmag_um_err = Kmag_err[unmatched_IR_ind]

#now plot
fig = plt.figure(figsize=(8,8))

print(len(Kmag_m), len(Kmag_um))

print(Kmag_um)
plt.errorbar(IR_m['ap_flux_B3'], Kmag_m, xerr=IR_m['ap_flux_err_B3'], yerr=Kmag_m_err, marker='o', linestyle='')
plt.errorbar(B3_upper_lim, Kmag_um, yerr=Kmag_um_err, marker='>', linestyle='')

plt.ylabel('K band magnitude')
plt.xlabel('Band 3 flux (mJy)')

plt.semilogx()

plt.savefig('/home/jotter/nrao/plots/Kband_B3_flux_plot.png', dpi=400)
plt.close()


MLLA_h_ind = np.where(np.isnan(IR_tab['Hmag1'])==False)
HC2000_h_ind = np.where(np.isnan(IR_tab['Hmag'])==False)

Hmag = np.repeat(np.nan, len(IR_tab))
Hmag_err = np.repeat(np.nan, len(IR_tab))

Hmag[HC2000_h_ind] = IR_tab['Hmag'][HC2000_h_ind]
Hmag_err[HC2000_h_ind] = IR_tab['e_Hmag'][HC2000_h_ind]
Hmag[MLLA_h_ind] = IR_tab['Hmag1'][MLLA_h_ind]
Hmag_err[MLLA_h_ind] = IR_tab['e_Hmag1'][MLLA_h_ind]

HK_color = Hmag-Kmag
HK_color_err = np.sqrt(Hmag_err**2 + Kmag_err**2)

HK_color_m = HK_color[matched_IR_ind]
HK_color_m_err = HK_color_err[matched_IR_ind]

HK_color_um = HK_color[unmatched_IR_ind]
HK_color_um_err = HK_color_err[unmatched_IR_ind]


fig = plt.figure(figsize=(8,8))

plt.errorbar(IR_m['ap_flux_B3'], HK_color_m, xerr=IR_m['ap_flux_err_B3'], yerr=HK_color_m_err, marker='o', linestyle='')
plt.errorbar(B3_upper_lim, HK_color_um, yerr=HK_color_um_err, marker='>', linestyle='')

plt.ylabel('H-K magnitude')
plt.xlabel('Band 3 flux (mJy)')

plt.semilogx()

plt.savefig('/home/jotter/nrao/plots/H-K_B3_flux_plot.png', dpi=400)
plt.close()
