from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from mpl_plot_templates import asinh_norm
from matplotlib.patches import Circle
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

IR = Table.read('/home/jotter/nrao/tables/A11_MLLA_r0.5_matched.fits')
data = Table.read('/home/jotter/nrao/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')

#compare non matches to our FOV
B3fl = fits.open('/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
img = B3fl[0].data.squeeze()
header = B3fl[0].header
wcs = WCS(header).celestial

TL_sky = wcs.all_pix2world(0,0, 0)
BR_sky = wcs.all_pix2world(len(img),len(img[0]), 0)

RA_upper = TL_sky[0]
RA_lower = BR_sky[0]
DEC_left = TL_sky[1]
DEC_right = BR_sky[1]

RA_ind1 = np.where(IR['RA'] <= RA_upper)[0]
RA_ind2 = np.where(IR['RA'] >= RA_lower)[0]

DEC_ind1 = np.where(IR['DEC'] <= DEC_right)[0]
DEC_ind2 = np.where(IR['DEC'] >= DEC_left)[0]

ra_ind = np.intersect1d(RA_ind1, RA_ind2)
dec_ind = np.intersect1d(DEC_ind1, DEC_ind2)
pos_ind = np.intersect1d(ra_ind, dec_ind)
IR = IR[pos_ind]

IR.write('/home/jotter/nrao/tables/A11_MLLA_r0.5_matched_FOV.fits')

nm_ind = np.where(np.isnan(IR['RA_B3']) == True)[0]
IR_nm = IR[nm_ind]

#next plot with circles

fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.15,0.1,0.8,0.8],projection=wcs)
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
dec.set_major_formatter('dd:mm:ss.s')
ra.ticklabels.set_fontsize(18)
dec.ticklabels.set_fontsize(18)

plt.imshow(img, origin='lower', transform=ax.get_transform(wcs), norm=asinh_norm.AsinhNorm(), vmin=-0.001, vmax=0.01)

IR_nm_pix = wcs.all_world2pix(IR_nm['RA']*u.degree, IR_nm['DEC']*u.degree, 0)
B3_pix = wcs.all_world2pix(data['RA_B3']*u.degree, data['DEC_B3']*u.degree, 0)

for ind in range(len(B3_pix[0])):
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=100, fill=False, color='red')
    ax.add_patch(circ)
    #ax.text(B3_pix[0][ind]-1, B3_pix[1][ind]+3, ind, color='red')

print(len(IR_nm))
    
for ind in range(len(IR_nm_pix[0])):
    circ = Circle((IR_nm_pix[0][ind], IR_nm_pix[1][ind]), radius=110, fill=False, color='orange')
    ax.add_patch(circ)
    #ax.text(IR_nm_pix[0][ind]-1, IR_nm_pix[1][ind]+3, f"{IR_nm['MLLA'][ind]},{IR_nm['Seq'][ind]}", color='blue')

ax.axis([0,len(img),0,len(img)])
plt.savefig(f'/home/jotter/nrao/plots/missed_IR_B3_overlay_tc.png',dpi=400)
plt.close()
