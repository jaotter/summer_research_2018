from fit_background import *
from astropy.io import fits
from astropy.table import Table

catalog = Table(fits.getdata('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_conv_snr.fits'))

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

B3_srcs = [2,17,21,29,46,49,57,6,7,13,19,26,8,69,70,71,72,73,74,75]
B3_xstddevs = [20,30,30,30,30,30,30,30,40,20,40,40,50,30,30,30,30,30,30,30]
B3_ystddevs = [20,30,30,30,30,30,30,30,40,20,40,40,50,30,30,30,30,30,30,30]
B3_zooms = [2,2,2,2,2,2,2,2,1.5,2,2,2,2,2,2,2,2.5,2.5,2.5,2]
B3_xmeans = [30,0,0,0,0,0,0,0,50,-50,-50,0,0,0,0,0,0,0,0,0]
B3_ymeans = [-30,0,0,0,0,0,0,0,0,-20,0,0,0,0,0,0,0,0,0,0]

B6_srcs = [2,3,7,8,16,17,19,24,34,40]
B6_xstddevs = [50,50,30,50,50,40,50,40,30,40]
B6_ystddevs = [50,50,30,50,50,40,50,40,30,60]
B6_zooms = [2,1.5,2,1.5,1.5,2,1.5,2,1.5,2]
B6_xmeans = [0,0,0,0,0,0,0,0,10,0]
B6_ymeans = [0,0,0,0,0,0,0,0,10,0]

B7_srcs = [0,3,16,17,27]
B7_xstddevs = [30,50,40,20,30]
B7_ystddevs = [30,50,40,20,30]
B7_xmeans = [-40,0,-25,20,0]
B7_ymeans = [40,0,25,-15,0]
B7_zooms = [1.3,1.5,2,3.5,2]

srcs = [B3_srcs, B6_srcs, B7_srcs]
xstddevs = [B3_xstddevs, B6_xstddevs, B7_xstddevs]
ystddevs = [B3_ystddevs, B6_ystddevs, B7_ystddevs]
zooms = [B3_zooms, B6_zooms, B7_zooms]
xmeans = [B3_xmeans, B6_xmeans, B7_xmeans]
ymeans = [B3_ymeans, B6_ymeans, B7_ymeans]

catalog['success'] = np.array(np.repeat('-', len(catalog)), dtype='S8')

mask_size = 1.5
for i in range(len(bands)): #loop thru bands
    for j in range(len(srcs[i])):
        print(srcs[i][j])

        if srcs[i][j] == 3:
            mask_size = 10
        fit = fit_source(srcs[i][j], imgs[i], names[i], bands[i], xstddevs[i][j], ystddevs[i][j], xmeans[i][j], ymeans[i][j], zoom=zooms[i][j], nonconv_img = nonconv_imgs[i], mask_size = mask_size)
        print(fit)
        cat_ind = np.where(catalog['D_ID'] == fit['D_ID'][0])[0]
        for nm in fit.colnames:
            catalog[nm][cat_ind] = fit[nm][0]
        mask_size = 1.5
catalog['pa_B3'] = catalog['pa_B3']%360-90
catalog['pa_B6'] = catalog['pa_B6']%360-90
catalog['pa_B7'] = catalog['pa_B7']%360-90

#ind1 = np.where(catalog['D_ID'] == 22)[0][0]
#ind2 = np.where(catalog['D_ID'] == 31)[0][0]
#ind3 = np.where(catalog['D_ID'] == 37)[0][0]
#ind4 = np.where(catalog['D_ID'] == 48)[0][0]
#catalog.remove_rows([ind1,ind2,ind3,ind4])

catalog.write('../tables/r0.5_catalog_conv_bgfitted_snr.fits',overwrite=True)

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
