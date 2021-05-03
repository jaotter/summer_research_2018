from fit_background import *
from astropy.io import fits
from astropy.table import Table


B3_img = '/lustre/cv/observers/cv-12578/orion_disks/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
B3_name = 'B3_huge_bg_cat'

B6_img = '/home/jotter/nrao/images/B6_convolved_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'
B6_name = 'B6_huge_bg_cat'
B7_img = '/home/jotter/nrao/images/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
B7_name = 'B7_conv_bg_cat'
B6_img_nonconv = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'
B7_img_nonconv = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'

fit_param_tab = Table.read('/lustre/cv/observers/cv-12578/orion_disks/summer_research_2018/tables/b3_fit_params_may21.csv')
#fit_param_tab = Table.read('/home/jotter/nrao/tables/b6_fit_params.csv')


#catalog = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_b6_catalog_may21.fits')
catalog = Table.read('/lustre/cv/observers/cv-12578/orion_disks/tables/r0.5_b3_catalog_may21.fits')

#catalog.rename_column('Seq_B3', 'Seq')
catalog.rename_column('D_ID', 'Seq')


for fit_row in fit_param_tab:
    #print(fit_row['source'])
    src = fit_row['\ufeffsource']
    print(src)
    zoom = float(fit_row['zoom'])
    if fit_row['notes'] == 'nonconv':
        img = B6_img_nonconv
    else:
        img = B3_img
    if fit_row['xmean'] == '-':
        bgfit = False
        xmean, ymean, xsigma, ysigma = [None,None,None,None]
    else:
        bgfit = True
        xmean, ymean, xsigma, ysigma = (float(fit_row['xmean']), float(fit_row['ymean']), float(fit_row['xsigma']), float(fit_row['ysigma']))
    fit = fit_source(src, img, B3_name, 'B3', fit_bg=bgfit, bg_stddev_x=xsigma, bg_stddev_y=ysigma, bg_mean_x=xmean, bg_mean_y=ymean, zoom=zoom, max_offset_in_beams=5)
    print(fit)
    cat_ind = np.where(catalog['Seq'] == fit['Seq'][0])[0]
    for nm in fit.colnames:
        catalog[nm][cat_ind] = fit[nm][0]

catalog['pa_B3'] = catalog['pa_B3']%360-90
#catalog['pa_B6'] = catalog['pa_B6']%360-90
#catalog = catalog[np.where(catalog['Seq'] == 117)]
print(catalog)

catalog.rename_column('Seq', 'B3_Seq')
#catalog.write('/home/jotter/nrao/summer_research_2018/tables/r0.5_b6_catalog_bgfit_may21.fits',overwrite=True)
catalog.write('/lustre/cv/observers/cv-12578/orion_disks/tables/r0.5_b3_catalog_bgfit_may21.fits',overwrite=True)
