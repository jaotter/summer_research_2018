from fit_background import *
from astropy.io import fits
from astropy.table import Table

catalog = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_jun20.fits')

bands = ['B3','B6', 'B7']
B6_img = '/home/jotter/nrao/images/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B6_name = 'B6_conv_bg_cat'
B7_img = '/home/jotter/nrao/images/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
B7_name = 'B7_conv_bg_cat'
B3_img = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
B3_name = 'B3_bg_r0.5.clean0.05mJy.allbaselines_cat'

B6_nonconv_img = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'
B7_nonconv_img = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'

nonconv_imgs = [None, B6_nonconv_img, B7_nonconv_img]

imgs = [B3_img, B6_img, B7_img]
names = [B3_name, B6_name, B7_name]

#dict order: source:[zoom, xmean, ymean, xsigma, ysigma]
B3_dict = {81:[5.5,0,0,20,20], 72:[5,0,0,30,30], 49:[2,0,0,30,30], 16:[2,0,0,30,30], 67:[2,0,0,30,30], 61:[2,0,0,30,30], 10:[2,0,0,30,30], 27:[4.5,0,0,30,30],
           29:[4,0,0,30,30], 34:[2,-50,-20,20,20], 39:[2,-50,0,40,40], 13:[2,0,0,40,40], 28:[2,0,0,50,50], 7:[2.5,50,50,0,0], 66:[2,0,0,30,30], 47:[2,0,0,30,30],
           21:[4,0,0,30,30], 73:[5.5,0,0,30,30], 4:[2.5,0,0,30,30], 2:[2,0,0,30,30], 3:[4,0,0,50,50], 20:[3,0,0,30,30], 37:[4.5,0,0,30,30], 44:[4.5,0,0,30,30], 47:[3,0,0,30,30],
           52:[4,0,0,30,30], 57:[4,0,0,30,30], 62:[4,0,0,30,30], 63:[5.5,0,0,30,30], 73:[5.5,0,0,30,30], 74:[5.5,0,0,30,30], 76:[7,0,0,30,30], 78:[4.5,0,0,30,30], 79:[4.5,0,0,30,30],
           80:[5.5,0,0,30,30], 83:[8,0,0,30,30], 84:[6,0,0,30,30]}

B6_dict = {81:[2,0,0,50,50], 29:[2,0,0,30,30], 28:[1.5,0,0,50,50], 72:[3.5,0,0,35,35], 39:[1.5,0,0,50,50], 18:[2,0,0,40,40], 50:[2,0,0,40,60], 10:[4.5,0,0,30,30], 13:[4,0,0,40,40],
           16:[12.7,0,0,30,30], 18:[4,0,0,30,30], 33:[4.5,0,0,30,30], 34:[8,0,0,30,30], 37:[3.5,0,0,30,30], 55:[3.5,0,0,30,30], 71:[9.5,0,0,30,30], 76:[3.5,0,0,30,30], 80:[8,0,0,30,30],
           83:[4,20,-20,30,30], 84:[4.5,0,0,30,30]}

B6_nonconv_img = [16,34,71,80,83] #B6 sources which should use non convolved image

B7_dict = {12:[1.3,-40,40,30,30], 37:[5,0,0,30,30], 14:[2,0,0,30,30], 16:[14,0,0,30,30], 18:[6,0,0,30,30], 27:[6,0,0,30,30], 29:[4,0,0,30,30], 33:[4,0,0,30,30], 34:[8,0,0,30,30],
           49:[3,0,0,30,30], 50:[8,0,0,30,30], 71:[5.5,0,0,30,30], 72:[5.5,0,0,30,30], 75:[4,0,0,30,30], 76:[8,0,0,30,30], 79:[3,0,0,30,30], 80:[12,0,0,30,30], 81:[12,0,0,30,30],
           83:[12,0,0,30,30]}

B7_nonconv_img = [16,18,33,34,50,71,76,80,81,83]


#catalog['success'] = np.array(np.repeat('-', len(catalog)), dtype='S8')

mask_size = 1.5

src_dict_list = [B3_dict, B6_dict, B7_dict]
nonconv_img_list = [[], B6_nonconv_img, B7_nonconv_img]

for i in range(len(src_dict_list)): #loop thru bands
    src_dict = src_dict_list[i]
    nonconv_sources = nonconv_img_list[i]
    
    for src in src_dict.keys():
        print(src) 
        
        if src == 3:
            mask_size = 10

        zoom, xmean, ymean, xsigma, ysigma = src_dict[src]

        image = imgs[i]
        if src in nonconv_sources:
            image = nonconv_imgs[i]
        
        fit = fit_source(src, image, names[i], bands[i], xsigma, ysigma, xmean, ymean, zoom=zoom, nonconv_img = nonconv_imgs[i], mask_size = mask_size, max_offset_in_beams=10)
        print(fit)
        cat_ind = np.where(catalog['D_ID'] == fit['D_ID'][0])[0]
        for nm in fit.colnames:
            catalog[nm][cat_ind] = fit[nm][0]
        mask_size = 1.5
        
catalog['pa_B3'] = catalog['pa_B3']%360-90
catalog['pa_B6'] = catalog['pa_B6']%360-90
catalog['pa_B7'] = catalog['pa_B7']%360-90

catalog.write('../tables/r0.5_catalog_bgfit_jun20.fits',overwrite=True)

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
