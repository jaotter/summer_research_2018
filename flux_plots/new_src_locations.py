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



tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim.fits')

b3fl = fits.open('/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits')
b3header = b3fl[0].header
b3data = b3fl[0].data
b3fl.close()

mywcs = WCS(b3header).celestial

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection=mywcs)
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
dec.set_major_formatter('dd:mm:ss.s')
ra.ticklabels.set_fontsize(18)
dec.ticklabels.set_fontsize(18)

ax.imshow(b3data, origin='lower', transform=ax.get_transform(mywcs)) #norm=asinh_norm.AsinhNorm()


new_srcs = np.array([8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124])

B3_pix = mywcs.all_world2pix(B3_table['RA_B3'][new_srcs]*u.degree, B3_table['DEC_B3'][new_srcs]*u.degree, 0)

#text_col = 'whitesmoke'


for ind, src in new_srcs:
    did = B3_table['ID'][src]
    col = 'tab:green'
    
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=10, fill=False, color=col)
    ax.add_patch(circ)
    ax.text(B3_pix[0][ind]-4, B3_pix[1][ind]+12, str(did), color='tab:green')

   
plt.savefig(f'/home/jotter/nrao/plots/new_src_loc.pdf',dpi=300,bbox_inches='tight')
plt.close()







