import numpy as np
import matplotlib.pyplot as plt
import radio_beam
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

FWHM_TO_SIGMA = 1/np.sqrt(8*np.log(2))

data = Table(fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_fixed.fits'))

band='B3'
#band = 'B6'

img_path = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
#img_path = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'

fl = fits.open(img_path)
header = fl[0].header
mywcs = WCS(header).celestial
beam = radio_beam.Beam.from_fits_header(header)
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg

deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_'+band]) == False)[0]

sigma_maj_deconv_pix = ((data['fwhm_maj_deconv_'+band]*u.arcsec).to(u.deg)/pixel_scale).decompose()*FWHM_TO_SIGMA
sigma_min_deconv_pix = ((data['fwhm_min_deconv_'+band]*u.arcsec).to(u.deg)/pixel_scale).decompose()*FWHM_TO_SIGMA

int_flux = 2*np.pi*data['gauss_amp_'+band]*sigma_maj_deconv_pix*sigma_min_deconv_pix
#int_flux = 2*np.pi*data['gauss_amp_'+band]*data['fwhm_maj_deconv_'+band]*data['fwhm_min_deconv_'+band]*FWHM_TO_SIGMA**2*1000
#int_flux2 = 2*np.pi*data['gauss_amp_'+band]*(data['fwhm_maj_deconv_'+band]*u.arcsec).to(u.degree)*(data['fwhm_min_deconv_'+band]*u.arcsec).to(u.degree)*FWHM_TO_SIGMA**2

ap_flux = data['ap_flux_'+band]
#gauss_flux = 2*np.pi*

#plt.plot([np.min(int_flux[deconv_ind]), np.max(int_flux[deconv_ind])], [np.min(int_flux[deconv_ind]), np.max(int_flux[deconv_ind])], color='k')
plt.scatter(int_flux[deconv_ind], ap_flux[deconv_ind])
plt.xlabel('integrated flux '+band+' (Jy)')
plt.ylabel('aperture flux ' +band+' (Jy)')
plt.ylim(-0.005, 0.04)
#plt.xlim(-0.005, 0.2)
plt.savefig('plots/'+band+'_int_ap_flux_comp.png')


