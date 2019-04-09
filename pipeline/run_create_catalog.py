from create_catalog import *


B3dr = '/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r'
#B3_imgs = [B3dr+'-2.clean0.5mJy.500klplus.image.tt0.pbcor.fits', B3dr+'0.5.clean1mJy.allbaselines.image.tt0.pbcor.fits']
#B3_img_names = ['B3_r-2.clean0.5mJy.500klplus', 'B3_r0.5.clean1mJy.allbaselines']

B6dr = '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r'
#B6_imgs = [B6dr+'-2.clean0.5mJy.500klplus.image.tt0.pbcor.fits', B6dr+'0.5.clean1mJy.allbaselines.image.tt0.pbcor.fits']
#B6_img_names = ['B6_r-2.clean0.5mJy.500klplus', 'B6_r0.5.clean1mJy.allbaselines']

B7dr = '/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r'
#B7_imgs = [B7dr+'-2.clean0.5mJy.500klplus.image.tt0.pbcor.fits', B7dr+'0.5.clean1mJy.allbaselines.image.tt0.pbcor.fits']
#B7_img_names = ['B7_r-2.clean0.5mJy.500klplus', 'B7_r0.5.clean1mJy.allbaselines']

B3_imgs = B3dr+'0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
B3_img_names = 'B3_r0.5.clean0.05mJy.allbaselines.deepmask'

#B6_imgs = B6dr+'0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B6_img_names = 'B6_r0.5.clean0.05mJy.150mplus.deepmask'

#B7_imgs = B7dr+'0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
B7_img_names = 'B7_r0.5.clean0.05mJy.250klplus.deepmask'

#single_img_catalog(B3_imgs, B3_img_names, B6_imgs, B6_img_names, B7_imgs, B7_img_names, 'r0.5_catalog_altB6')

##line to create catalog with convolved data
B3_img = B3_imgs
B3_name = B3_img_names

B6_img = '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
#B6_img = B6dr+'0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B6_name = 'conv_'+B6_img_names
#B6_name = B6_img_names

B7_img = '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
#B7_img = B7dr+'0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
B7_name = 'conv_'+B7_img_names
#B7_name = B7_img_names

single_img_catalog(B3_img, B3_name, B6_img, B6_name, B7_img, B7_name, 'r0.5_catalog_convB6B7_ppbeamtest')
