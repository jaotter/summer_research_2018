from fit_background import *

src = 20
fit_bg = False

bg_xmean = 0
bg_ymean = 0
bg_xsigma = 3
bg_ysigma = 3
zoom = 2

band = 'B3'

B3_img = '/lustre/cv/observers/cv-12578/orion_disks/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
B3_name = 'B3_huge_bg'


#B6_img = '/home/jotter/nrao/images/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#B6_name = 'B6_conv_r0.5.clean.0.05mJy.150mplus'
#B7_img = '/home/jotter/nrao/images/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
#B7_name = 'B7_conv_r0.5.clean.0.05mJy.250klplus'

#B6_img = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#B6_name = 'B6_bg_r0.5.clean.0.05mJy.150mplus'
#B7_img = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
#B7_name = 'B7_bg_r0.5.clean.0.05mJy.250klplus'

srcs = [23,25,29,33,36,40,42,43,45,52,54,55,56,57,58,59,61,63,64,65,66,68,69,70,71,74,76,77,79,80,81,83,84,85,86,87,88,89,92,93,94,95,97,98,100,102,103,105,106,107,108,109,110,111,112,117,119,120,121,122,123,124,125]

for src in srcs:
    fit = fit_source(src, B3_img, B3_name, band, fit_bg=fit_bg, bg_stddev_x=bg_xsigma, bg_stddev_y=bg_ysigma, bg_mean_x=bg_xmean, bg_mean_y=bg_ymean, zoom=3, max_offset_in_beams=5)


#fit = fit_source(src, B3_img, B3_name, band, fit_bg=fit_bg, bg_stddev_x=bg_xsigma, bg_stddev_y=bg_ysigma, bg_mean_x=bg_xmean, bg_mean_y=bg_ymean, zoom=zoom, max_offset_in_beams=5)
#print(fit['fwhm_maj_deconv_B6'], fit['fwhm_maj_deconv_err_B6'])
#print(fit['RA_B6'], fit['DEC_B6'])
print(fit)
#print(fit['ap_flux_B6'])
#fit params: - default xmean 0, ymean 0, zoom 1
#B7:
#src 0, x 31, y 31, zoom 1.3, xmean -40, ymean 40
#src 3, x 50, y 50, zoom 1.5
#src 16, x 40, y 40, zoom 2, xmean -25, ymean 25
#src 17, x 40, y 40, zoom 1.5 - still kinda bad

#B6:
#src 2, x 50, y 50, zoom 2 
#src 3, x 50, y 50, zoom 1.5
#src 7, x 30, y 30, zoom 2
#src 8, x 50, y 50, zoom 1.5
#src 16, x 50, y 50, zoom 1.5
#src 17, x 40, y 40, zoom 2
#src 19, x 50, y 50, zoom 1.5
#src 24, x 40, y 40, zoom 2
#src 34, x 31, y 31, zoom 1.5, xmean 10, ymean 10 - still kinda bad
#src 40, x 40, y 60, zoom 2
