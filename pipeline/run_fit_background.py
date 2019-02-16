from fit_background import *

bands = ['B6', 'B7']
B6_img = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B6_name = 'B6_conv_r0.5.clean.0.05mJy.150mplus'
B7_img = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
B7_name = 'B7_conv_r0.5.clean.0.05mJy.250klplus'

imgs = [B6_img, B7_img]
names = [B6_name, B7_name]

B6_srcs = [2,3,7,8,16,17,19,24,34,40]
B6_xstddevs = [50,50,30,50,50,40,50,40,30,40]
B6_ystddevs = [50,50,30,50,50,40,50,40,30,60]
B6_zooms = []

for i in range(len(bands)):
    


fit = fit_source(13, B6_img, B6_name, band, 50, 50, 0, 0, zoom=1.5)
print(fit)

#fit params: - default xmean 0, ymean 0, zoom 1
#B7:
#src 0, x 31, y 31, zoom 1.3, xmean -40, ymean 40
#src 3, x 50, y 50, zoom 1.5
#src 16, x 40, y 40, zoom 2, xmean -25, ymean 25
#src 17 - need more than zoom in?

#B6:
#src 2, x 50, y 50, zoom 2 
#src 3, x 50, y 50, zoom 1.5
#src 7, x 30, y 30, zoom 2
#src 8, x 50, y 50, zoom 1.5
#src 16, x 50, y 50, zoom 1.5
#src 17, x 40, y 40, zoom 2
#src 19, x 50, y 50, zoom 1.5
#src 24, x 40, y 40, zoom 2
#src 34, x 31, y 31, zoom 1.5, xmean 10, ymean 10
#src 40, x 40, y 60, zoom 2
