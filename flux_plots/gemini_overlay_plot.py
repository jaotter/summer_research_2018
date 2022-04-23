from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
from matplotlib.patches import Circle, RegularPolygon
from mpl_plot_templates import asinh_norm
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np

plt.ion()

basepath = '/Users/adam/work/students/JustinOtter/summer_research_2018/'

img_path = f'{basepath}/Trapezium_GEMS_mosaic_redblueorange_normed_small_contrast_bright_photoshop.png'
im = img.imread(img_path)

B3_table = Table.read(f'{basepath}/tables/r0.5_catalog_bgfit_may21_ulim.fits')

header = fits.Header.fromtextfile(f'{basepath}/trapezium_small_photoshop.wcs')
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

b7_sources = np.where(B3_table['SNR_B7'] > 3)[0]
b6_sources = np.where(B3_table['SNR_B6'] > 3)[0]

b6_sources = np.setdiff1d(b6_sources, b7_sources)

b3_sources = np.setdiff1d(np.arange(len(B3_table)),b6_sources)
b3_sources = np.setdiff1d(b3_sources,b7_sources)

print("lengths: ", len(b7_sources), len(b6_sources), len(b3_sources))


B3_table = Table.read(f'{basepath}/final_tables/datafile4.txt', format='ascii.cds')
Fbid = B3_table['ID'][~B3_table['F2016'].mask]
coupid = B3_table['ID'][~B3_table['COUP'].mask]
IRid = B3_table['ID'][~B3_table['MLLA'].mask]
new_srcs = B3_table['ID'][B3_table['f_ID'] == 'a']
print(f"new sources: {len(new_srcs)}, coup: {len(coupid)}, forbrich: {len(Fbid)}, IR: {len(IRid)}")

IRsel = np.isin(B3_table['ID'], IRid)
XRsel = np.isin(B3_table['ID'], coupid)
radsel = np.isin(B3_table['ID'], Fbid) & ~IRsel
newsel = np.isin(B3_table['ID'], new_srcs)
mmsel = (~IRsel) & (~XRsel) & (~radsel) & (~newsel)

coords = SkyCoord((B3_table['RAs'].quantity + (5*u.hour).to(u.s) + (35*u.min).to(u.s))*(15*u.arcsec/u.s),
                  -5*u.deg - B3_table['DEm'].quantity - B3_table['DEs'].quantity, frame='icrs', )
IRscat = ax.scatter(coords.ra[IRsel], coords.dec[IRsel], transform=ax.get_transform('icrs'), facecolors='none', edgecolor='lime', marker='p')
XRscat = ax.scatter(coords.ra[XRsel], coords.dec[XRsel], transform=ax.get_transform('icrs'), facecolors='none', edgecolor='cyan', marker='d')
radscat = ax.scatter(coords.ra[radsel], coords.dec[radsel], transform=ax.get_transform('icrs'), facecolors='none', edgecolor='orange', marker='^')
mmscat = ax.scatter(coords.ra[mmsel], coords.dec[mmsel], transform=ax.get_transform('icrs'), facecolors='none', edgecolor='blue', marker='o')


ax.set_ylim(0,1250)
ax.set_xlim(100,1556)

ax.set_xlabel('Right Ascension', fontsize=18)
ax.set_ylabel('Declination', fontsize=18)

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_ir_xray_radio_past.png',dpi=200,bbox_inches='tight')


ax.scatter(coords.ra[newsel], coords.dec[newsel], transform=ax.get_transform('icrs'), facecolors='none', c='red', marker='*')

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_ir_xray_radio_past_withNew.png',dpi=200,bbox_inches='tight')


B3_coord = SkyCoord(ra=83.8104625, dec=-5.37515556, unit=u.degree)
B3_coord_pix = mywcs.all_world2pix(B3_coord.ra, B3_coord.dec, 0)
B3_radius = 137.6*u.arcsecond/2
B3_coord_R = SkyCoord(ra=B3_coord.ra, dec=B3_coord.dec+B3_radius)
B3_R_pix = mywcs.all_world2pix(B3_coord_R.ra, B3_coord_R.dec,0)
B3_radius_pix = B3_coord_pix[1] - B3_R_pix[1]

B3_circ = Circle((B3_coord_pix[0], B3_coord_pix[1]), radius=B3_radius_pix, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='--')
ax.add_patch(B3_circ)

B6_coord = SkyCoord(ra=83.8104625, dec=-5.37515556, unit=u.degree)
B6_coord_pix = mywcs.all_world2pix(B6_coord.ra, B6_coord.dec, 0)
B6_radius = 28.6*u.arcsecond
B6_coord_R = SkyCoord(ra=B6_coord.ra, dec=B6_coord.dec+B6_radius)
B6_R_pix = mywcs.all_world2pix(B6_coord_R.ra, B6_coord_R.dec,0)
B6_radius_pix = B6_coord_pix[1] - B6_R_pix[1]

B6_circ = Circle((B6_coord_pix[0], B6_coord_pix[1]), radius=B6_radius_pix, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='--')
ax.add_patch(B6_circ)

B7_coord = SkyCoord(83.8104626, -5.37515542, unit=u.degree)
B7_coord_pix = mywcs.all_world2pix(B7_coord.ra, B7_coord.dec, 0)
B7_radius = 12.7*u.arcsec
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

ax.set_ylim(0,1250)
ax.set_xlim(100,1556)
#ax.set_ylim(BL_pix[1], TR_pix[1])
#ax.set_xlim(BL_pix[0], TR_pix[0])

#plt.tight_layout()

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_wavedet.png',dpi=200,bbox_inches='tight')


for scat in (IRscat, radscat, XRscat, mmscat):
    scat.set_visible(False)

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_newOnly.png',dpi=200,bbox_inches='tight')

plt.ion()

lores_cont = fits.open('/Users/adam/Dropbox/Orion_ALMA/FITS/Orion_NW_12m_continuum.fits')
ax.contour(lores_cont[0].data, transform=ax.get_transform(WCS(lores_cont[0].header)), levels=[0.1, 0.2, 0.5, 0.75, 1], colors=['lime']*6, linewidth=0.75)

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_newOnly_contours.png',dpi=200,bbox_inches='tight')

ax.axis((750.3124848440887, 928.9302949291337, 508.5377590382199, 694.7765956868934))

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_newOnly_contours_zoom.png',dpi=200,bbox_inches='tight')

for scat in (IRscat, radscat, XRscat, mmscat):
    scat.set_visible(True)

plt.savefig(f'{basepath}/plots/gemini_B3_overlay_contours_zoom.png',dpi=200,bbox_inches='tight')
