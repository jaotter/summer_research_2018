from astropy.io import fits
from astropy.table import Table
from radio_beam import Beam

import numpy as np
import astropy.units as u

B3img = fits.open('/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
B3beam = Beam.from_fits_header(B3img[0].header)

data = Table.read('../tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')

nondeconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B3']) == True)[0]

fwhm_UL = []
for ind in nondeconv_ind:
    src_size = Beam(major=(data['fwhm_maj_B3'][ind]+data['fwhm_maj_err_B3'][ind])*u.arcsec, minor=(data['fwhm_min_B3'][ind]+data['fwhm_min_err_B3'][ind])*u.arcsec, pa=(data['pa_B3'][ind]-90)*u.degree)
    try:
        deconv_UL = src_size.deconvolve(B3beam)
        fwhm_UL.append(deconv_UL.major.value)
    except ValueError:
        fwhm_UL.append(np.nan)
        print('source '+str(data['D_ID'][ind])+' could not be deconvolved')

print(fwhm_UL)
