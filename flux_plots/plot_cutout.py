from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from mpl_plot_templates import asinh_norm
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord


def plot_cutout(img_path, coord, size, savefile):
    fl = fits.open(img_path)
    img = fl[0].data
    wcs = WCS(fl[0].header).celestial
    cutout = Cutout2D(img, coord, size, wcs)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_axes([0.15,0.1,0.8,0.8],projection=wcs)
    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm:ss.s')
    ra.ticklabels.set_fontsize(18)
    dec.ticklabels.set_fontsize(18)

    plt.imshow(cutout.data, transform=ax.get_transform(wcs), norm=asinh_norm.AsinhNorm())
    plt.colorbar()
    
    plt.savefig(f'/home/jotter/nrao/plots/IR_cutouts/{savefile}.png',dpi=400)
    plt.close()


img = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
IR = Table.read('/home/jotter/nrao/tables/A11_MLLA_r0.5_HC2000_LRY2000_matched_FOV.fits')

um_ind = np.where(np.isnan(IR['RA_B3']) == True)

IR_um = IR[um_ind]

IR_um = IR_um[23::]

um_coords = SkyCoord(ra=IR_um['RAJ2000_1'], dec=IR_um['DEJ2000_1'], unit=u.degree)
print(IR_um['MLLA'].data)

for ind, coord in enumerate(um_coords):
    srcname = IR_um['MLLA'][ind]
    print(srcname)
    plot_cutout(img, coord, 2*u.arcsec, f'src_{srcname}MLLA_cutout.png')

#src1ind = np.where(IR['MLLA'] == 506)[0]
#src2ind = np.where(IR['MLLA'] == 628)[0]
#680, 638, 671, 614, 606, 603, 610, 591, 630, 639, 664
