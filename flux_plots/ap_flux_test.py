from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import numpy as np
import radio_beam

FWHM_TO_SIGMA = 

B6_nonconv_img = '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'

nonconv_fl = fits.open(B6_nonconv_img)
nonconv_header = nonconv_fl[0].header
nonconv_wcs = WCS(nonconv_header).celestial

nonconv_beam = Beam.from_fits_header(nonconv_header)
pixel_scale_nonconv = np.abs(nonconv_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
ppbeam_nonconv = (nonconv_beam.sr/(pixel_scale_nonconv**2)).decompose().value

data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_fixed.fits')

ind = np.where(data['D_ID'] == 11)[0]

center_coord = SkyCoord(data['RA_B6'][ind]*u.deg, data['DEC_B6']*u.deg, frame='icrs')
center_coord_pix = center_coord.to_pixel(nonconv_wcs)
center_coord_pix_reg = region.PixCoord(center_coord_pix[0], center_coord_pix[1])

pix_maj_fwhm_nonconv = ((data['fwhm_maj_deconv_B6'][ind]*u.arcsec).to(u.deg)/pixel_scale).decompose()
pix_min_fwhm_nonconv = ((data['fwhm_min_deconv_B6'][ind]*u.arcsec).to(u.deg)/pixel_scale).decompose()
