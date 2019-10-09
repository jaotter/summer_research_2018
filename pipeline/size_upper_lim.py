from astropy.io import fits
from radio_beam import Beam
from astropy.convolution import convolve_fft
from astropy.wcs import WCS
from matplotlib import patches
from astropy.coordinates import SkyCoord, Angle
from gaussfit_catalog import gaussfit_catalog
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import regions
import sys

B3file = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
B6file = ''
B7file = ''

orig_fl = fits.open(B3file)
orig_data = fits.getdata(B3file)

#create fake image data - uniform disk with radius r
r = 0.0005 #r is radius in arcseconds

fake_data = np.zeros(orig_data.size)
fake_mgrid = np.mgrid[:len(orig_data), :len(orig_data)]
print(fake_data.size, fake_mgrid.size)

sys.exit('done')

new_fl = orig_fl
new_fl[0].data = fake_img
wcs = WCS(new_fl[0].header).celestial
B3beam = Beam.from_fits_header(new_fl[0].header)

B3_pixel_scale = np.abs(wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
B3_kernel = B3beam.as_kernel(B3_pixel_scale)

convolved_img = convolve_fft(fake_img, B3_kernel)
new_fl[0].data = convolved_img
new_fl.writeto('/users/jotter/summer_research_2018/fake_conv_img.fits', overwrite=True)

loc = wcs.wcs_pix2world(500,500,1)
reg = regions.CircleSkyRegion(center=SkyCoord(loc[0], loc[1], unit='deg'), radius=0.5*u.arcsec, meta={'text':'test'})
reg_pix = reg.to_pixel(wcs)

gaussfit = gaussfit_catalog('/users/jotter/summer_research_2018/fake_conv_img.fits', [reg], Angle(0.2, 'arcsecond'), savepath='/users/jotter/summer_research_2018/')

source_size = Beam(major=gaussfit['test']['fwhm_major'], minor=gaussfit['test']['fwhm_minor'], pa=(gaussfit['test']['pa'].value+90)*u.degree)
try:
    deconv_size = source_size.deconvolve(B3beam)
    print('gaussfit deconv major: '+str(gaussfit['test']['deconv_fwjm_major'])+' minor: '+str(gaussfit['test']['deconv_fwhm_minor'])+' deconv + 90 degrees major: '+str(deconv_size.major.value)+' minor: '+str(deconv_size.minor.value))

except ValueError:
    print('could not be deconvolved')
    
'''
plt.imshow(convolved_img, cmap='viridis', origin='lower', interpolation='nearevvvvvst')

center_deg = SkyCoord(gaussfit['test']['center_x'], gaussfit['test']['center_y'], unit='deg')
center_pix = skycoord_to_pixel(center_deg, wcs)

fwhm_maj_deconv_pix = ((gaussfit['test']['deconv_fwhm_major']*u.arcsec).to(u.degree)/pixel_scale).decompose()
fwhm_min_deconv_pix = ((gaussfit['test']['deconv_fwhm_minor']*u.arcsec).to(u.degree)/pixel_scale).decompose()
ellipse = patches.Ellipse(center_pix, width=fwhm_min_deconv_pix, height=fwhm_maj_deconv_pix, angle=gaussfit['test']['deconv_pa'])

plt.add_patch(ellipse)

plt.savefig('/users/jotter/summer_research_2018/fake_img_fit.png')
'''
