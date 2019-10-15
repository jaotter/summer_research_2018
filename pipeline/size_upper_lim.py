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

def draw_circle(r,size,flux): #draw circle radius r in 2d rectangular array arr
    xx, yy = np.mgrid[:size, :size]
    radius = np.sqrt((xx - size/2)**2 + (yy - size/2)**2)
    circle_arr = np.zeros((size, size))
    circle_ind = np.where(radius < r)
    circle_arr[circle_ind] = flux
    return circle_arr

    
B3file = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
B6file = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B7file = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'

B3noise = 5e-5
B6noise = 1e-4
B7noise = 1e-3

B3flux = 1e-4
B6flux = 0
B7flux = 0

img = B3file
flux = B3flux #flux of each pixel in fake img
noise = B3noise #noise std dev

radii_au = [1,2,3,4,5,8,10,15,20,30,40,50]*u.AU
radii_as = (radii_au.to(u.pc)/(400*u.pc)*u.rad).to(u.arcsecond).value
radii_au = radii_au.value
#radii_as = [0.2,0.15,0.1,0.08,0.06,0.05,0.045,0.04,0.035,0.03,0.025]

orig_fl = fits.open(img)

wcs = WCS(orig_fl[0].header).celestial
beam = Beam.from_fits_header(orig_fl[0].header)

pixel_scale = np.abs(wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg

kernel = beam.as_kernel(pixel_scale)

size=1000 #square size of fake image

recovered_fwhm_maj = []

for i,rad in enumerate(radii_as):
  #create fake image data - uniform disk with radius r
  r_arcsec = rad*u.arcsec #radius in arcsec
  r = (r_arcsec/pixel_scale).decompose() #radius in pixel
  fake_data = draw_circle(r,size,flux)

  convolved_fake_img = convolve_fft(fake_data, kernel)

  #add noise similar to images
  noise_arr = np.random.normal(0, noise, (size,size))

  noisy_conv_fake_img = noise_arr + convolved_fake_img

  new_fl = orig_fl
  new_fl[0].data = noisy_conv_fake_img
  new_fl.writeto('/users/jotter/summer_research_2018/pipeline/size_lims/fake_conv_img_r%s.fits'%(str(radii_au[i])), overwrite=True)


  loc = wcs.wcs_pix2world(size/2,size/2,1)
  reg = regions.CircleSkyRegion(center=SkyCoord(loc[0], loc[1], unit='deg'), radius=0.5*u.arcsec, meta={'text':'r=%s'%str(radii_au[i])})
  reg_pix = reg.to_pixel(wcs)

  gaussfit = gaussfit_catalog('/users/jotter/summer_research_2018/pipeline/size_lims/fake_conv_img_r%s.fits'%(str(radii_au[i])), [reg], Angle(0.3, 'arcsecond'), savepath='/users/jotter/summer_research_2018/pipeline/size_lims/', max_radius_in_beams=15)

  source_size = Beam(major=gaussfit['r=%s'%str(radii_au[i])]['fwhm_major'], minor=gaussfit['r=%s'%str(radii_au[i])]['fwhm_minor'], pa=(gaussfit['r=%s'%str(radii_au[i])]['pa'].value+90)*u.degree)
  try:
      deconv_size = source_size.deconvolve(beam)
      print(deconv_size)
      recovered_fwhm_maj.append(deconv_size.major.value)
  except ValueError:
      print('could not be deconvolved')
      recovered_fwhm_maj.append(0)



recovered_fwhm_maj_au = (((recovered_fwhm_maj*u.arcsecond).to(u.radian)).value*(400*u.pc)).to(u.AU)
plt.figure()
plt.scatter(radii_au, recovered_fwhm_maj_au)
plt.xlabel('disk radius (AU)')
plt.ylabel('recovered major FWHM (AU)')
#plt.xlim(0,50)
#plt.ylim(0,50)
plt.plot([0,50],[0,50])

plt.savefig('size_lims/recov_radius.png')
