from fit_background import *

band = 'B6'
B6_img = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B6_name = 'B6_conv_r0.5.clean.0.05mJy.150mplus'
B7_img = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
B7_name = 'B7_conv_r0.5.clean.0.05mJy.250klplus'

fit = fit_source(16, B7_img, B7_name, band, 30, 30, 50, 0, zoom=1.25)
print(fit)

#fit params:
#B7:
#src 0, x 30, y30 - also try zoom in
#src 16 - need zoom in
#src 17 - need more than zoom in?

#B6:
#src 2, x 50, y 30
#src 7 - not working
#src 16, x 80, y 80
#src 17, x 50, y 100
#src 19, x 50, y 50
#src 24, x 70, y 40
#src 34 - not working - try to zoom in
#src 40, x 40, y 60
