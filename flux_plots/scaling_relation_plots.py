from astropy.table import Table
from astropy.io import fits

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

band='B3'

#three plots: R and L, R and M, L and M
data = Table.read('../tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)

Rarcsec = data['fwhm_maj_deconv_'+band]*u.arcsec
d = (414*u.pc).to(u.au)
Rau = Rarcsec.to(u.rad)*d
print(Rau)

