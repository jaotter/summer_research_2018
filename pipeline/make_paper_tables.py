from astropy.table import Table
from astropy.io import fits

import radio_beam.Beam
import numpy as np
import astropy.units as u

data = Table.read('../tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')

#calculate quantities for each band
bands = ['B3', 'B6', 'B7']
imgs = ['/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
nonconv_imgs = ['/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']

for b in range(len(bands)):
    band = bands[b]
    fl = fits.open(imgs[b])



table1 = Table(data['D_ID'], data['RA_B3'], data['DEC_B3'], data['RA_err_B3'], data['DEC_err_B3'],
               data['ap_flux_B3'], data['ap_flux_err_B3'], data['ap_flux_B6'], data['ap_flux_err_B6'],
               data['ap_flux_B7'], data['ap_flux_err_B7'], data['gauss_amp_B3'],
               data['gauss_amp_err_B3'], data['gauss_amp_B6'], data['gauss_amp_err_B6'],
               data['gauss_amp_B7'], data['gauss_amp_err_B7'], data['fwhm_maj_deconv_B3'],
               data['fwhm_maj_deconv_err_B3'], data['fwhm_min_B3'],
               data['fwhm_min_err_B3'], data['pa_B3'], data['pa_err_B3'],
               data['fwhm_maj_B6'], data['fwhm_maj_err_B6'], data['fwhm_min_B6'],
               data['fwhm_min_err_B6'], data['pa_B6'], data['pa_err_B6'],
               data['fwhm_maj_B7'], data['fwhm_maj_err_B7'], data['fwhm_min_B7'],
               data['fwhm_min_err_B7'], data['pa_B7'], data['pa_err_B7'],
               data['fwhm_maj_deconv_B3'], data['fwhm_min_deconv_B3'], data['pa_deconv_B3'], 
               data['fwhm_maj_deconv_B6'], data['fwhm_min_deconv_B6'], data['pa_deconv_B6'], 
               data['fwhm_maj_deconv_B7'], data['fwhm_min_deconv_B7'], data['pa_deconv_B7'])


