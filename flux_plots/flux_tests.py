import numpy as np
import matplotlib.pyplot as plt
import radio_beam
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

FWHM_TO_SIGMA = 1/np.sqrt(8*np.log(2))

data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_nonconv_apflux_final.fits')

data_conv = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')



#band='B3'
band = 'B6'
#band = 'B7'

#img_path_conv = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
img_path_conv = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#img_path_conv = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'


#img_path = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
img_path = '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#img_path = '/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'


fl = fits.open(img_path)
header = fl[0].header
mywcs = WCS(header).celestial
beam = radio_beam.Beam.from_fits_header(header)
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
ppbeam = (beam.sr/((pixel_scale)**2)).decompose().value

flux_ind = np.where(np.isnan(data['ap_flux_'+band]) == False)[0]

sigma_maj_deconv_pix = ((data['fwhm_maj_'+band]*u.arcsec).to(u.deg)/pixel_scale).decompose()*FWHM_TO_SIGMA
sigma_min_deconv_pix = ((data['fwhm_min_'+band]*u.arcsec).to(u.deg)/pixel_scale).decompose()*FWHM_TO_SIGMA

int_flux = 2*np.pi*data['gauss_amp_'+band]*sigma_maj_deconv_pix*sigma_min_deconv_pix/ppbeam

fl_conv = fits.open(img_path_conv)
header_conv = fl_conv[0].header
mywcs_conv = WCS(header_conv).celestial
beam_conv = radio_beam.Beam.from_fits_header(header_conv)
pixel_scale_conv = np.abs(mywcs_conv.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
ppbeam_conv = (beam_conv.sr/((pixel_scale_conv)**2)).decompose().value

flux_ind_conv = np.where(np.isnan(data_conv['ap_flux_'+band]) == False)[0]

sigma_maj_deconv_pix_conv = ((data_conv['fwhm_maj_'+band]*u.arcsec).to(u.deg)/pixel_scale_conv).decompose()*FWHM_TO_SIGMA
sigma_min_deconv_pix_conv = ((data_conv['fwhm_min_'+band]*u.arcsec).to(u.deg)/pixel_scale_conv).decompose()*FWHM_TO_SIGMA

#int_flux_conv = 2*np.pi*data_conv['gauss_amp_'+band]*sigma_maj_deconv_pix_conv*sigma_min_deconv_pix_conv/ppbeam_conv
int_flux_conv = 2*np.pi*data_conv['gauss_amp_'+band]*sigma_maj_deconv_pix_conv*sigma_min_deconv_pix_conv/ppbeam

#int_flux = ((2*np.pi*data['gauss_amp_'+band]*data['fwhm_maj_deconv_'+band]*u.arcsec*data['fwhm_min_deconv_'+band]*u.arcsec*(FWHM_TO_SIGMA**2))/beam.sr).decompose()

ap_flux = data['ap_flux_'+band]
ap_flux_conv = data_conv['ap_flux_'+band]

#plt.plot([np.min(int_flux[deconv_ind]), np.max(int_flux[deconv_ind])], [np.min(int_flux[deconv_ind]), np.max(int_flux[deconv_ind])], color='k')

#print(ap_flux[11], int_flux[11])
#plt.scatter(int_flux[deconv_ind], ap_flux[deconv_ind])

full_ind = np.intersect1d(flux_ind, flux_ind_conv)

flux_ind_conv = np.delete(flux_ind_conv, np.where(flux_ind_conv==3)[0])
flux_ind_conv = np.delete(flux_ind_conv, np.where(flux_ind_conv==13)[0])

#plt.scatter(int_flux[full_ind], int_flux_conv[full_ind])
plt.scatter(int_flux_conv[flux_ind_conv], ap_flux_conv[flux_ind_conv])
#plt.scatter(int_flux[deconv_ind], ap_flux[deconv_ind])

plt.plot([0,1], [0,1], color='k')
plt.xlabel('integrated flux (convolved) '+band+' (Jy)')
#plt.ylabel('integrated flux (convolved) '+band+' (Jy)')
plt.ylabel('aperture flux (convolved)' +band+' (Jy)')
plt.ylim(-0.01, 0.1)
plt.xlim(-0.01, 0.1)
plt.savefig('plots/'+band+'_int_ap_flux_comp_conv.png')


