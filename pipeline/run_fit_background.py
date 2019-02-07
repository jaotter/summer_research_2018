from fit_background import *

band = 'B6'
B6_img = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B6_name = 'B6_conv_r0.5,clean.0.05mJy.150mplus'
#B7_img = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
#B7_name = 'conv_maxroffset2_'+B7_img_names

fit = fit_source(B6_img, B6_name, band, 50, 50)
print(fit)
