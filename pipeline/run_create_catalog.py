from create_catalog import *


B3_img = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
B3_name = 'B3_r0.5.clean0.05mJy.allbaselines.deepmask'

#B6_img = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#B6_name = 'B6_r0.5.clean0.05mJy.150mplus.deepmask'

#B7_img = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
#B7_name = 'B7_r0.5.clean0.05mJy.250klplus.deepmask'

#B6_img_conv = '/home/jotter/nrao/images/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#B6_name_conv = 'conv_B6_r0.5.clean0.05mJy.150mplus.deepmask'

#B7_img_conv = '/home/jotter/nrao/images/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
#B7_name_conv = 'conv_B7_r0.5.clean0.05mJy.250klplus.deepmask'

#single_img_catalog(B3_img, B3_name, B6_img_conv, B6_name_conv, B7_img_conv, B7_name_conv, 'r0.5_catalog_feb21', nonconv_B6_img=B6_img, nonconv_B7_img=B7_img)

#single_img_catalog(B3_img, B3_name, B6_img, B6_name, B7_img, B7_name, 'r0.5_catalog_nonconv_apflux_final')


huge_B3_catalog('/lustre/cv/observers/cv-12578/orion_disks/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits', 'r0.5_b3_catalog_may21')
