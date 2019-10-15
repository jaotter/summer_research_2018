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

def draw_circle(r,size,total_flux): #draw circle radius r in 2d rectangular array arr
    xx, yy = np.mgrid[:size, :size]
    radius = np.sqrt((xx - size/2)**2 + (yy - size/2)**2)
    circle_arr = np.zeros((size, size))
    circle_ind = np.where(radius < r)
    circle_arr[circle_ind] = total_flux
    return circle_arr

    
B3file = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
B6file = ''
B7file = ''

orig_fl = fits.open(B3file)

wcs = WCS(orig_fl[0].header).celestial
beam = Beam.from_fits_header(orig_fl[0].header)

pixel_scale = np.abs(wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg

#create fake image data - uniform disk with radius r
r_arcsec = 0.1*u.arcsec #radius in arcsec
r = (r_arcsec/pixel_scale).decompose() #radius in pixel
size = 1000 #size of square fake data in pixels
tot_flux = 1e-4 #total flux of source
fake_data = draw_circle(r,size,tot_flux)

plt.imshow(fake_data)
plt.savefig('size_lims/r0.1ascircle.png')

kernel = beam.as_kernel(pixel_scale)
convolved_fake_img = convolve_fft(fake_data, kernel)

#next step: add noise similar to images
noise_rms = 5e-5 
noise_arr = np.random.normal(0, noise_rms, (size,size))

noisy_conv_fake_img = noise_arr + convolved_fake_img

new_fl = orig_fl
new_fl[0].data = noisy_conv_fake_img
new_fl.writeto('/users/jotter/summer_research_2018/pipeline/size_lims/fake_conv_img.fits', overwrite=True)

sys.exit('done')

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
