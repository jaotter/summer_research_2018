from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
from matplotlib.patches import Circle
from mpl_plot_templates import asinh_norm
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


cube_path = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
E18_table_path = '/home/jotter/nrao/tables/eis_coord_table.fits'
E18_table = Table.read(E18_table_path)

B3_table = Table.read('/home/jotter/nrao/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')

fl = fits.open(cube_path)
img = fl[0].data.squeeze()
header = fl[0].header
mywcs = WCS(header).celestial

fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.15,0.1,0.8,0.8],projection=mywcs)
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
dec.set_major_formatter('dd:mm:ss.s')
ra.ticklabels.set_fontsize(18)
dec.ticklabels.set_fontsize(18)

plt.imshow(img, origin='lower', transform=ax.get_transform(mywcs), norm=asinh_norm.AsinhNorm(), vmin=-0.001, vmax=0.01)

E18_pix = mywcs.all_world2pix(E18_table['RA']*u.degree, E18_table['DEC']*u.degree, 0)
B3_pix = mywcs.all_world2pix(B3_table['RA_B3']*u.degree, B3_table['DEC_B3']*u.degree, 0)

print(len(img), len(img[0]))
print(B3_pix[0], B3_pix[1])

for ind in range(len(B3_pix[0])):
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=100, fill=False, color='red')
    ax.add_patch(circ)
    #ax.text(table['pixel_x'][i]-1, table['pixel_y'][i]+3, src_ID, color='red')

for ind in range(len(E18_pix[0])):
    circ = Circle((E18_pix[0][ind], E18_pix[1][ind]), radius=110, fill=False, color='orange')
    ax.add_patch(circ)

ax.axis([0,len(img),0,len(img)])
plt.savefig(f'/home/jotter/nrao/plots/E18_B3_overlay.png',dpi=400)
plt.close()
