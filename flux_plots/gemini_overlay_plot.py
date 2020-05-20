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

B3_table = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_apr20.fits')

fl = fits.open(cube_path)
img = fl[0].data.squeeze()
header = fl[0].header
mywcs = WCS(header).celestial

TR_pix = mywcs.all_world2pix(83.805*u.degree, -5.3715*u.degree, 0)
cut_img = img[:int(TR_pix[1]), :int(TR_pix[0])]

fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.15,0.1,0.8,0.8],projection=mywcs)
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
#ra.set_major_formatter('d.ddd')
dec = ax.coords['dec']
dec.set_major_formatter('dd:mm:ss.s')
#dec.set_major_formatter('d.ddd')
ra.ticklabels.set_fontsize(18)
dec.ticklabels.set_fontsize(18)

ax.imshow(cut_img, origin='lower', transform=ax.get_transform(mywcs), norm=asinh_norm.AsinhNorm(), vmin=-0.001, vmax=0.005)

E18_pix = mywcs.all_world2pix(E18_table['RA']*u.degree, E18_table['DEC']*u.degree, 0)
B3_pix = mywcs.all_world2pix(B3_table['RA_B3']*u.degree, B3_table['DEC_B3']*u.degree, 0)

for ind in range(len(B3_pix[0])):
    if B3_pix[0][ind] < TR_pix[0] and B3_pix[1][ind] < TR_pix[1]:
        circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=100, fill=False, color='red')
        ax.add_patch(circ)
        #ax.text(B3_pix[0][ind]-1, B3_pix[1][ind]+3, B3_table['D_ID'][ind], color='red')

for ind in range(len(E18_pix[0])):
    circ = Circle((E18_pix[0][ind], E18_pix[1][ind]), radius=110, fill=False, color='orange')
    ax.add_patch(circ)

#ax.axis([0,len(img),0,len(img)])
plt.savefig(f'/home/jotter/nrao/plots/E18_comp/E18_B3_overlay.png',dpi=300)
plt.close()



