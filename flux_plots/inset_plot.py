from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.gridspec import GridSpec
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from radio_beam import Beam
from matplotlib.patches import Ellipse, Rectangle
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

def inset_plot_3band(page_num):
    #create inset plot for sources detected in 3 bands
    tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
    b3path = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
    b6path = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'
    b7path = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'

    b3_fl = fits.open(b3path)
    b3_data = b3_fl[0].data
    b3_wcs = WCS(b3_fl[0].header).celestial
    b3beam = Beam.from_fits_header(b3_fl[0].header)
    b3_fl.close()
    b6_fl = fits.open(b6path)
    b6_data = b6_fl[0].data
    b6_wcs = WCS(b6_fl[0].header).celestial
    b6beam = Beam.from_fits_header(b6_fl[0].header)
    b6_fl.close()
    b7_fl = fits.open(b7path)
    b7_data = b7_fl[0].data
    b7_wcs = WCS(b7_fl[0].header).celestial
    b7beam = Beam.from_fits_header(b6_fl[0].header)
    b7_fl.close()

    wcs_list = [b3_wcs, b6_wcs, b7_wcs]
    data_list = [b3_data, b6_data, b7_data]
    beam_list = [b3beam, b6beam, b7beam]
    
    srcind = np.where(np.isnan(tab['ap_flux_B7']) == False)[0] #36 sources detected in B7
    nondet = [14, 37, 38]
    nondet_bands = np.array([[1], [1,2], [1,2]], dtype='object')
    srcind = np.sort(np.concatenate((srcind, [37,38])))
    
    if page_num == 1:
        srcind = srcind[:12]
    if page_num == 2:
        srcind = srcind[12:24]
    if page_num == 3:
        srcind = srcind[24:]

    tab = tab[srcind]

    new_srcs = [8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124]
    
    fig = plt.figure(figsize=(8,9))
    outer_gs = GridSpec(6, 2, figure=fig, wspace=0.05, hspace=0.01)

    ind = 0
    for gs_row in range(6):
        for gs_col in range(2):

            #now create inner gridspec for 3 bands
            inner_gs = outer_gs[gs_row, gs_col].subgridspec(1,3,wspace=0)
            wavelength = ['3mm', '1.3mm', '0.85mm']
            if ind >= len(srcind):
                continue
            
            for i in range(3):
                img = data_list[i]
                mywcs = wcs_list[i]
                ax = fig.add_subplot(inner_gs[0,i], projection=mywcs)
                src_coord = SkyCoord(ra=tab['RA_B3'][ind], dec=tab['DEC_B3'][ind], unit=u.degree)
                coord_pix = src_coord.to_pixel(mywcs)
                pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
                size = (0.5*u.arcsecond).to(u.degree)/pixel_scale
                size = size.decompose().value
                cutout = Cutout2D(img, coord_pix, size, mywcs, mode='partial')
                ax.imshow(cutout.data, origin='lower', cmap='viridis', transform=ax.get_transform(mywcs))

                if i == 0:
                    bbox_params = {'visible':False}
                    if tab['ID'][ind] in new_srcs:
                        bbox_params = {'edgecolor':'red','fill':None,'boxstyle':'round'}
                    ax.text(10,100,tab['ID'][ind],color='white', bbox=bbox_params)
                    ax.tick_params(axis='x', which='both', direction='in', color='black', labelbottom=False)
                    ax.tick_params(axis='y', which='both', direction='in', color='black', labelleft=False)

                else:
                    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labeltop=False)
                    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)

                ax.set_xlabel('')
                ax.set_ylabel('')

                if tab['ID'][ind] in nondet:
                    nondet_band_inds = nondet_bands[np.where(tab['ID'][ind] == nondet)[0]]
                    for band_ind in nondet_band_inds[0]:
                        if i == band_ind:
                            ax.text(70,10,'S/N<3',color='white', fontsize=8)
                    
                if gs_row == 0 and gs_col == 0:
                    spacings = [70,60,50]
                    ax.text(spacings[i],100,wavelength[i],color='white', fontsize=8)
                    ax.tick_params(axis='both', which='both', direction='in', color='black')

                    beam = beam_list[i]
                    ellipse_center = mywcs.pixel_to_world(20,20)
                    r = Ellipse((ellipse_center.ra.value, ellipse_center.dec.value),
                            beam.major/(1. * u.deg), beam.minor/(1. * u.deg), float(beam.pa/(1. * u.deg)), edgecolor='yellow',
                            facecolor='yellow', alpha=0.5, transform=ax.get_transform('icrs'))
                    ax.add_patch(r)

                    scalebar_deg = ((50/400) * u.arcsecond).to(u.degree)
                    
                    if i == 2:
                        rect_cent = mywcs.pixel_to_world(115,5)
                        r = Rectangle((rect_cent.ra.value, rect_cent.dec.value), scalebar_deg.value, scalebar_deg.value/5,
                                      edgecolor='black', facecolor='gray', transform=ax.get_transform('icrs'))
                        ax.add_patch(r)
                        ax.text(75, 15, '50 AU', fontsize=8)
    
                                        
            ind += 1

    plt.savefig(f'/home/jotter/nrao/plots/insets/b7_page{page_num}.pdf', bbox_inches='tight')


def inset_plot_2band(page_num):
    #create inset plot for sources detected in 3 bands
    tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
    b3path = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
    b6path = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'
    
    b3_fl = fits.open(b3path)
    b3_data = b3_fl[0].data
    b3_wcs = WCS(b3_fl[0].header).celestial
    b3beam = Beam.from_fits_header(b3_fl[0].header)
    b3_fl.close()
    b6_fl = fits.open(b6path)
    b6_data = b6_fl[0].data
    b6_wcs = WCS(b6_fl[0].header).celestial
    b6beam = Beam.from_fits_header(b6_fl[0].header)
    b6_fl.close()
    
    wcs_list = [b3_wcs, b6_wcs]
    data_list = [b3_data, b6_data]
    beam_list = [b3beam, b6beam]
    
    b7_srcind = np.where(np.isnan(tab['ap_flux_B7']) == False)[0] #34 sources detected in B7
    b6_srcind = np.where(np.isnan(tab['ap_flux_B6']) == False)[0]
    srcind = np.setdiff1d(b6_srcind, b7_srcind) #16 sources in b6 and not b7
    nondet = [3,19,43,53,58]
    srcind = np.sort(np.concatenate((srcind, nondet))) #22 sources total
    print(len(srcind))

    tab = tab[srcind]

    new_srcs = [8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124]
    
    fig = plt.figure(figsize=(7.5,10.5))
    outer_gs = GridSpec(8, 3, figure=fig, wspace=0.05, hspace=0.01)

    ind = 0
    for gs_row in range(8):
        for gs_col in range(3):

            #now create inner gridspec for 3 bands
            inner_gs = outer_gs[gs_row, gs_col].subgridspec(1,2,wspace=0)
            wavelength = ['3mm', '1.3mm']
            if ind >= len(srcind):
                continue
            
            for i in range(2):
                img = data_list[i]
                mywcs = wcs_list[i]
                ax = fig.add_subplot(inner_gs[0,i], projection=mywcs)
                src_coord = SkyCoord(ra=tab['RA_B3'][ind], dec=tab['DEC_B3'][ind], unit=u.degree)
                coord_pix = src_coord.to_pixel(mywcs)
                pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
                size = (0.5*u.arcsecond).to(u.degree)/pixel_scale
                size = size.decompose().value
                cutout = Cutout2D(img, coord_pix, size, mywcs, mode='partial')
                ax.imshow(cutout.data, origin='lower', cmap='viridis', transform=ax.get_transform(mywcs))
                if i == 0:
                    bbox_params = {'visible':False}
                    if tab['ID'][ind] in new_srcs:
                        bbox_params = {'edgecolor':'red','fill':None,'boxstyle':'round'}
                    ax.text(10,100,tab['ID'][ind],color='white', bbox=bbox_params)
                    ax.tick_params(axis='x', which='both', direction='in', color='black', labelbottom=False)
                    ax.tick_params(axis='y', which='both', direction='in', color='black', labelleft=False)

                else:
                    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labeltop=False)
                    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)

                ax.set_xlabel('')
                ax.set_ylabel('')


                if tab['ID'][ind] in nondet:
                    if i == 1:
                        ax.text(70,10,'S/N<3',color='white', fontsize=8)
                
                
                if gs_row == 0 and gs_col == 0:
                    spacings = [70,60,50]
                    ax.text(spacings[i],100,wavelength[i],color='white', fontsize=8)
                    ax.tick_params(axis='both', which='both', direction='in', color='black')

                    beam = beam_list[i]
                    ellipse_center = mywcs.pixel_to_world(20,20)
                    r = Ellipse((ellipse_center.ra.value, ellipse_center.dec.value),
                            beam.major/(1. * u.deg), beam.minor/(1. * u.deg), float(beam.pa/(1. * u.deg)), edgecolor='yellow',
                            facecolor='yellow', alpha=0.5, transform=ax.get_transform('icrs'))
                    ax.add_patch(r)

                    scalebar_deg = ((50/400) * u.arcsecond).to(u.degree)
                    #scalebar_1kpc_pix = float(scalebar_100AU_arcsec/(0.5*u.arcsec) * 1.*u.kpc)
                    
                    if i == 1:
                        rect_cent = mywcs.pixel_to_world(115,25)
                        r = Rectangle((rect_cent.ra.value, rect_cent.dec.value), scalebar_deg.value, scalebar_deg.value/5,
                                      edgecolor='black', facecolor='gray', transform=ax.get_transform('icrs'))
                        ax.add_patch(r)
                        ax.text(75, 35, '50 AU', fontsize=8)
    
                        #highlight newly added sources?
                    #add axes labels
                    
                    
            ind += 1

    plt.savefig(f'/home/jotter/nrao/plots/insets/b6_page{page_num}.pdf', bbox_inches='tight')


def inset_plot_1band(page_num):
    #create inset plot for sources detected in 3 bands
    tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
    b3path = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
    b6path = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'
    
    b3_fl = fits.open(b3path)
    b3_data = b3_fl[0].data
    b3_wcs = WCS(b3_fl[0].header).celestial
    b3beam = Beam.from_fits_header(b3_fl[0].header)
    b3_fl.close()
    
    b7_srcind = np.where(np.isnan(tab['ap_flux_B7']) == False)[0] #34 sources detected in B7
    b6_srcind = np.where(np.isnan(tab['ap_flux_B6']) == False)[0]
    b3b6_srcind = np.setdiff1d(np.arange(len(tab)), b6_srcind)
    srcind = np.setdiff1d(b3b6_srcind, b7_srcind)
    prev_nondet = [3,19,37,38,43,53,58]
    srcind = np.setdiff1d(srcind, prev_nondet)
    
    print(len(srcind))
    
    if page_num == 1:
        srcind = srcind[:36]
    if page_num == 2:
        srcind = srcind[36:]

        
    tab = tab[srcind]

    new_srcs = [8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124]
    
    fig = plt.figure(figsize=(8,9))
    outer_gs = GridSpec(7, 6, figure=fig, wspace=0.05, hspace=0.05)

    ind = 0
    for gs_row in range(6):
        for gs_col in range(6):

            #now create inner gridspec for 3 bands
            wavelength = '3mm'
            if ind >= len(srcind):
                continue
            
            img = b3_data
            mywcs = b3_wcs
            ax = fig.add_subplot(outer_gs[gs_row,gs_col], projection=mywcs)
            src_coord = SkyCoord(ra=tab['RA_B3'][ind], dec=tab['DEC_B3'][ind], unit=u.degree)
            coord_pix = src_coord.to_pixel(mywcs)
            pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
            size = (0.5*u.arcsecond).to(u.degree)/pixel_scale
            size = size.decompose().value
            cutout = Cutout2D(img, coord_pix, size, mywcs, mode='partial')
            ax.imshow(cutout.data, origin='lower', cmap='viridis', transform=ax.get_transform(mywcs))
            bbox_params = {'visible':False}
            if tab['ID'][ind] in new_srcs:
                bbox_params = {'edgecolor':'red','fill':None,'boxstyle':'round'}
            ax.text(10,100,tab['ID'][ind],color='white', bbox=bbox_params)
            ax.tick_params(axis='x', which='both', direction='in', color='black', labelbottom=False)
            ax.tick_params(axis='y', which='both', direction='in', color='black', labelleft=False)

            ax.set_xlabel('')
            ax.set_ylabel('')


            if gs_row == 0 and gs_col == 0:
                spacings = 70
                ax.text(spacings,100,wavelength,color='white', fontsize=8)
                ax.tick_params(axis='both', which='both', direction='in', color='black')

                beam = b3beam
                ellipse_center = mywcs.pixel_to_world(20,20)
                r = Ellipse((ellipse_center.ra.value, ellipse_center.dec.value),
                        beam.major/(1. * u.deg), beam.minor/(1. * u.deg), float(beam.pa/(1. * u.deg)), edgecolor='yellow',
                        facecolor='yellow', alpha=0.5, transform=ax.get_transform('icrs'))
                ax.add_patch(r)

                scalebar_deg = ((50/400) * u.arcsecond).to(u.degree)
                #scalebar_1kpc_pix = float(scalebar_100AU_arcsec/(0.5*u.arcsec) * 1.*u.kpc)

                
                rect_cent = mywcs.pixel_to_world(115,5)
                r = Rectangle((rect_cent.ra.value, rect_cent.dec.value), scalebar_deg.value, scalebar_deg.value/5,
                              edgecolor='black', facecolor='gray', transform=ax.get_transform('icrs'))
                ax.add_patch(r)
                ax.text(75, 15, '50 AU', fontsize=8)

                    
            ind += 1

    plt.savefig(f'/home/jotter/nrao/plots/insets/b3_page{page_num}.pdf', bbox_inches='tight')



#inset_plot_1band(1)
#inset_plot_1band(2)

inset_plot_2band(1)
    
#inset_plot_3band(1)
#inset_plot_3band(2)
#inset_plot_3band(3)
