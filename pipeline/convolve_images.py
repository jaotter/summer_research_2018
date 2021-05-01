import radio_beam

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import astropy.convolution as convolution
import astropy.units as u

B3file = fits.open('/lustre/cv/observers/cv-12578/orion_disks/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits')
#B3file = fits.open('/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits')

B3header = B3file[0].header
#B3img = B3file[0].data
B3beam = radio_beam.Beam.from_fits_header(B3header)

B6file = fits.open('/lustre/cv/observers/cv-12578/orion_disks/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits')
#B6file = fits.open('/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits')

B6header = B6file[0].header
B6img = B6file[0].data
B6beam = radio_beam.Beam.from_fits_header(B6header)
B6WCS = WCS(B6header).celestial
B6_pixel_scale = np.abs(B6WCS.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg

#B7file = fits.open('/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits')
#B7header = B7file[0].header
#B7img = B7file[0].data
#B7beam = radio_beam.Beam.from_fits_header(B7header)
#B7WCS = WCS(B7header).celestial
#B7_pixel_scale = np.abs(B7WCS.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg

B6_deconv = B3beam.deconvolve(B6beam)
#B7_deconv = B3beam.deconvolve(B7beam)

B6_kernel = B6_deconv.as_kernel(B6_pixel_scale)
#B7_kernel = B7_deconv.as_kernel(B7_pixel_scale)

convolved_B6 = convolution.convolve_fft(B6img, B6_kernel, allow_huge = True)
#convolved_B7 = convolution.convolve_fft(B7img, B7_kernel, allow_huge = True)

B6file[0].data = convolved_B6
B6file[0].header = B3beam.attach_to_header(B6file[0].header)


#B6file.writeto('/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', overwrite=True)
B6file.writeto('/lustre/cv/observers/cv-12578/orion_disks/B6_convolved_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits', overwrite=True)

#B7file[0].data = convolved_B7
#B7file[0].header = B3beam.attach_to_header(B7file[0].header)

#B7file.writeto('/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits', overwrite=True)
