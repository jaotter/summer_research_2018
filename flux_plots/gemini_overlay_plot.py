from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
from matplotlib.patches import Circle, Rectangle
from mpl_plot_templates import asinh_norm
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np


img_path = '/home/jotter/nrao/images/Trapezium_GEMS_mosaic_redblueorange_normed_small_contrast_bright_photoshop.png'
im = img.imread(img_path)

B3_table = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_feb21_ulim.fits')

header = fits.Header.fromtextfile('/home/jotter/nrao/images/trapezium_small_photoshop.wcs')
#header = fits.Header.fromtextfile('/home/jotter/nrao/images/fullimage.wcs')
mywcs = WCS(header).celestial

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection=mywcs)
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
#ra.set_major_formatter('d.ddd')
dec = ax.coords['dec']
dec.set_major_formatter('dd:mm:ss.s')
#dec.set_major_formatter('d.ddd')
ra.ticklabels.set_fontsize(18)
dec.ticklabels.set_fontsize(18)

ax.imshow(np.flip(np.rot90(im,2,axes=(0,1)),axis=1), origin='lower', transform=ax.get_transform(mywcs)) #norm=asinh_norm.AsinhNorm()

B3_pix = mywcs.all_world2pix(B3_table['RA_B3']*u.degree, B3_table['DEC_B3']*u.degree, 0)

inset1 = [30, 43, 31, 23, 38, 32, 28, 45, 20, 71, 22, 36]
inset2 = [34, 39, 41, 76, 79, 29, 80, 75, 35, 42, 49, 50] 
inset3 = [14, 81, 17, 16, 12, 11, 18, 27, 83, 33] 
insetb7 = np.concatenate((inset1, inset2, inset3))
inset4 = [10, 8, 13, 84, 54, 55, 53, 46]
inset5 = [5, 6, 9, 15, 40, 51, 21, 24, 48, 52, 63, 56, 61, 58, 7, 73, 4, 2, 3, 19, 44, 47, 78, 57, 60, 62, 59, 67, 66, 64, 70]

text_col = 'whitesmoke'

for ind in range(len(B3_pix[0])):
    did = B3_table['Seq'][ind]
    #if did in insetb7:
    #    col = 'tab:red'
    #if did in inset4:
    #    col = 'tab:pink'
    #if did in inset5:
    #    col = 'tab:green'
    col='tab:red'
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=10, fill=False, color=col)
    ax.add_patch(circ)
    #if did == 60 or did == 58 or did == 9 or did == 71 or did == 23 or did == 37 or did == 22 or did == 29 or did == 80 or did == 55:
    #    ax.text(B3_pix[0][ind]-5, B3_pix[1][ind]-18, str(did), color=text_col)
    #elif did == 45 or did == 27 or did == 64 or did == 70:
    #    ax.text(B3_pix[0][ind]-4, B3_pix[1][ind]+12, str(did), color='gray')
    #elif did == 43:
    #    ax.text(B3_pix[0][ind]-5, B3_pix[1][ind]-18, str(did), color='gray')
    #elif did == 20:
    #    ax.text(B3_pix[0][ind]-22, B3_pix[1][ind]-2, str(did), color='lightgray')
    #elif did == 50:
    #    ax.text(B3_pix[0][ind]-22, B3_pix[1][ind]-2, str(did), color=text_col)
    #elif did == 32 or did == 38:
    #    ax.text(B3_pix[0][ind]+10, B3_pix[1][ind], str(did), color=text_col)
    #else:
    #    ax.text(B3_pix[0][ind]-4, B3_pix[1][ind]+12, str(did), color=text_col)
    
    #ax.text(B3_pix[0][ind]-1, B3_pix[1][ind]+3, B3_table['D_ID'][ind], color='red')

B3_coord = SkyCoord(ra=83.8104625, dec=-5.37515556, unit=u.degree)
B3_coord_pix = mywcs.all_world2pix(B3_coord.ra, B3_coord.dec, 0)
B3_radius = 137.6*u.arcsecond/2
B3_coord_R = SkyCoord(ra=B3_coord.ra, dec=B3_coord.dec+B3_radius)
B3_R_pix = mywcs.all_world2pix(B3_coord_R.ra, B3_coord_R.dec,0)
B3_radius_pix = B3_coord_pix[1] - B3_R_pix[1]

B3_circ = Circle((B3_coord_pix[0], B3_coord_pix[1]), radius=B3_radius_pix, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='--')
ax.add_patch(B3_circ)

B6_coord = SkyCoord(ra=83.81446244, dec=-5.37913764, unit=u.degree)
B6_coord_pix = mywcs.all_world2pix(B6_coord.ra, B6_coord.dec, 0)
B6_length = 28.672*u.arcsec
B6_coord_BR = SkyCoord(ra=B6_coord.ra, dec=B6_coord.dec+B6_length)
B6_BR_pix = mywcs.all_world2pix(B6_coord_BR.ra, B6_coord_BR.dec,0)
B6_width = B6_BR_pix[1] - B6_coord_pix[1]
B6_coord_TL = SkyCoord(ra=B6_coord.ra+B6_length, dec=B6_coord.dec)
B6_TL_pix = mywcs.all_world2pix(B6_coord_TL.ra, B6_coord_TL.dec,0)
B6_height = B6_coord_pix[0] - B6_TL_pix[0]

B6_rect = Rectangle((B6_coord_pix[0], B6_coord_pix[1]), width=B6_width, height=B6_height, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='--')
ax.add_patch(B6_rect)

B7_coord = SkyCoord(83.8104626, -5.37515542, unit=u.degree)
B7_coord_pix = mywcs.all_world2pix(B7_coord.ra, B7_coord.dec, 0)
B7_radius = 12.9*u.arcsec
B7_coord_R = SkyCoord(ra=B7_coord.ra, dec=B7_coord.dec+B7_radius)
B7_R_pix = mywcs.all_world2pix(B7_coord_R.ra, B7_coord_R.dec,0)
B7_radius_pix = B7_coord_pix[1] - B7_R_pix[1]

B7_circ = Circle((B7_coord_pix[0], B7_coord_pix[1]), radius=B7_radius_pix, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='--')
ax.add_patch(B7_circ)

#ax.set_ylim(200,1000)
#ax.set_xlim(300,1200)

#BL_frame = SkyCoord(ra='5:35:17', dec='-5:23:10', unit=(u.hourangle,u.degree))
#BL_pix = mywcs.all_world2pix(BL_frame.ra.degree, BL_frame.dec.degree, 0)
#TR_frame = SkyCoord(ra='5:35:12', dec='-5:21:50', unit=(u.hourangle,u.degree))
#TR_pix = mywcs.all_world2pix(TR_frame.ra.degree, TR_frame.dec.degree, 0)

ax.set_ylim(0,1350)
ax.set_xlim(100,1556)
#ax.set_ylim(BL_pix[1], TR_pix[1])
#ax.set_xlim(BL_pix[0], TR_pix[0])

plt.tight_layout()
    
plt.savefig(f'/home/jotter/nrao/plots/gemini_B3_overlay.pdf',dpi=300,bbox_inches='tight')
plt.close()







