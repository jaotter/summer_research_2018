import numpy as np
import matplotlib.pyplot as plt
import radio_beam
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

FWHM_TO_SIGMA = 1/np.sqrt(8*np.log(2))

data_conv = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_fixed.fits')

data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_apflux_fixed.fits')

#band='B3'
band = 'B6'
#band = 'B7'

#img_path = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
img_path_conv = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#img_path_conv = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'

img_path = '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#img_path = '/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'


fl = fits.open(img_path)
header = fl[0].header
mywcs = WCS(header).celestial
beam = radio_beam.Beam.from_fits_header(header)
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
ppbeam = (beam.sr/((pixel_scale)**2)).decompose().value

fl_conv = fits.open(img_path_conv)
header_conv = fl[0].header
mywcs_conv = WCS(header_conv).celestial
beam_conv = radio_beam.Beam.from_fits_header(header_conv)
pixel_scale_conv = np.abs(mywcs_conv.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
ppbeam_conv = (beam_conv.sr/((pixel_scale_conv)**2)).decompose().value

deconv_ind = np.where(np.isnan(data_conv['fwhm_maj_deconv_'+band]) == False)[0]

#sigma_maj_deconv_pix = ((data['fwhm_maj_deconv_'+band]*u.arcsec).to(u.deg)/pixel_scale).decompose()*FWHM_TO_SIGMA
#sigma_min_deconv_pix = ((data['fwhm_min_deconv_'+band]*u.arcsec).to(u.deg)/pixel_scale).decompose()*FWHM_TO_SIGMA

#int_flux = 2*np.pi*data['gauss_amp_'+band]*sigma_maj_deconv_pix*sigma_min_deconv_pix/ppbeam

ap_flux = data['ap_flux_'+band]
ap_flux_conv = data_conv['ap_flux_'+band]

plt.scatter(ap_flux_conv[deconv_ind], ap_flux[deconv_ind])
plt.xlabel('aperture flux (convolved image) '+band+' (Jy)')
plt.ylabel('aperture flux ' +band+' (Jy)')
plt.ylim(-0.001, 0.03)
plt.xlim(-0.001, 0.03)
plt.savefig('plots/'+band+'_ap_flux_conv_comp.png')


