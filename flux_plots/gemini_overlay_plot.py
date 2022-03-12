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


img_path = '/home/jotter/nrao/images/Trapezium_GEMS_mosaic_redblueorange_normed_small_contrast_bright_photoshop.png'
im = img.imread(img_path)

B3_table = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim.fits')

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

b7_sources = np.where(B3_table['SNR_B7'] > 3)[0]
b6_sources = np.where(B3_table['SNR_B6'] > 3)[0]

b6_sources = np.setdiff1d(b6_sources, b7_sources)

b3_sources = np.setdiff1d(np.arange(len(B3_table)),b6_sources)
b3_sources = np.setdiff1d(b3_sources,b7_sources)

print(len(b7_sources), len(b6_sources), len(b3_sources))


text_col = 'whitesmoke'

irtab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
IRid = irtab['ID'].data

couptab = Table.read('/home/jotter/nrao/summer_research_2018/tables/COUP_r0.5_may21.fits')
coupid = couptab['ID']

Fbtab = Table.read('/home/jotter/nrao/summer_research_2018/tables/Forbrich2016_r0.5_may21.fits')
Fbid = Fbtab['ID']

new_srcs = [8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124]

for ind in range(len(B3_pix[0])):
    did = B3_table['ID'][ind]
    #if did in b7_sources:
    #    col = 'tab:red'
    #if did in b6_sources:
    #    col = 'tab:pink'
    #if did in b3_sources:
    #    col = 'tab:green'
    #col='tab:red'

    if did in IRid:
        #col = 'tab:pink'
        col = 'mistyrose'
        pentagon = RegularPolygon((B3_pix[0][ind], B3_pix[1][ind]), numVertices=5, radius=13, fill=False, color=col)
        pentagon.set_path_effects([PathEffects.Stroke(linewidth=2, foreground="black"), PathEffects.Normal()])
        ax.add_patch(pentagon)
        
    if did in coupid:
        #col = 'tab:green'
        col = 'honeydew'
        square = RegularPolygon((B3_pix[0][ind], B3_pix[1][ind]), numVertices=4, radius=8, fill=False, color=col)
        square.set_path_effects([PathEffects.Stroke(linewidth=2, foreground="black"), PathEffects.Normal()])
        ax.add_patch(square)
    if did in Fbid:
        #col = 'tab:orange'
        col = 'oldlace'
        triangle = RegularPolygon((B3_pix[0][ind], B3_pix[1][ind]), numVertices=3, radius=4, fill=False, color=col)
        triangle.set_path_effects([PathEffects.Stroke(linewidth=2, foreground="black"), PathEffects.Normal()])
        ax.add_patch(triangle)
    if did in new_srcs:
        #col = 'tab:red'
        col = 'azure'
        #circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=10, fill=False, color=col)    
        #ax.add_patch(circ)
        #ax.plot(B3_pix[0][ind], B3_pix[1][ind], marker = 'x', linestyle='')

        cent_x = B3_pix[0][ind]
        cent_y = B3_pix[1][ind]
        Path = mpath.Path
        path_data = [(Path.MOVETO, [cent_x-5, cent_y+5]),
                     (Path.LINETO, [cent_x+5, cent_y-5]),
                     (Path.MOVETO, [cent_x-5, cent_y-5]),
                     (Path.LINETO, [cent_x+5, cent_y+5])]
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path, fill=False, color=col)
        patch.set_path_effects([PathEffects.Stroke(linewidth=2, foreground="black"), PathEffects.Normal()])
        ax.add_patch(patch)

        
    if did not in IRid and did not in coupid and did not in Fbid and did not in new_srcs:
        #col = 'tab:blue'
        col = 'lavenderblush'
        circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=8, fill=False, color=col)
        circ.set_path_effects([PathEffects.Stroke(linewidth=2, foreground="black"), PathEffects.Normal()])
        ax.add_patch(circ)
        
        
    #circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=10, fill=False, color=col)    
    #ax.add_patch(circ)

        
    #else:
    #    col = 'tab:green'
    
    
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

ax.set_xlabel('RA', fontsize=18)
ax.set_ylabel('Declination', fontsize=18)
plt.tight_layout()
    
plt.savefig(f'/home/jotter/nrao/plots/gemini_B3_overlay_wavedet.pdf',dpi=300,bbox_inches='tight')
plt.close()







