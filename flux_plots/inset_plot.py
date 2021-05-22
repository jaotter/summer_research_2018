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
    
    srcind = np.where(np.isnan(tab['ap_flux_B7']) == False)[0] #34 sources detected in B7
    if page_num == 1:
        srcind = srcind[:12]
    if page_num == 2:
        srcind = srcind[12:24]
    if page_num == 3:
        srcind = srcind[24:]
        
    tab = tab[srcind]

    fig = plt.figure(figsize=(8,9))
    outer_gs = GridSpec(6, 2, figure=fig, wspace=0.05, hspace=0.01)

    ind = 0
    for gs_row in range(6):
        for gs_col in range(2):

            #now create inner gridspec for 3 bands
            inner_gs = outer_gs[gs_row, gs_col].subgridspec(1,3,wspace=0)
            #axs = plt.add_subplot(outer_gs[gs_row, gs_col])
            #axs = inner_gs.subplots()#subplot_kw={'projection':wcs_list[0]})
            #ax0 = fig.add_subplot(outer_gs[gs_row, gs_col])

            wavelength = ['3mm', '1.3mm', '0.85mm']
            
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
                    ax.text(10,100,tab['ID'][ind],color='white')
                    ax.tick_params(axis='x', which='both', direction='in', color='black', labelbottom=False)
                    ax.tick_params(axis='y', which='both', direction='in', color='black', labelleft=False)

                else:
                    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labeltop=False)
                    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)

                ax.set_xlabel('')
                ax.set_ylabel('')

                
                if gs_row == 0 and gs_col == 0:
                    spacings = [70,60,50]
                    ax.text(spacings[i],100,wavelength[i],color='white')
                    ax.tick_params(axis='both', which='both', direction='in', color='black')

                    beam = beam_list[i]
                    ellipse_center = mywcs.pixel_to_world(20,20)
                    r = Ellipse((ellipse_center.ra.value, ellipse_center.dec.value),
                            beam.major/(1. * u.deg), beam.minor/(1. * u.deg), float(beam.pa/(1. * u.deg)), edgecolor='yellow',
                            facecolor='yellow', alpha=0.5, transform=ax.get_transform('icrs'))
                    ax.add_patch(r)

                    scalebar_deg = ((50/400) * u.arcsecond).to(u.degree)
                    #scalebar_1kpc_pix = float(scalebar_100AU_arcsec/(0.5*u.arcsec) * 1.*u.kpc)
                    
                    if i == 2:
                        rect_cent = mywcs.pixel_to_world(115,5)
                        r = Rectangle((rect_cent.ra.value, rect_cent.dec.value), scalebar_deg.value, scalebar_deg.value/5,
                                      edgecolor='black', facecolor='gray', transform=ax.get_transform('icrs'))
                        ax.add_patch(r)
                        ax.text(75, 15, '50 AU', fontsize=8)
    

                    
                    #add beam
                    #add axes labels
                    
                    
            ind += 1

    plt.savefig(f'/home/jotter/nrao/plots/insets/b7_page{page_num}.pdf', bbox_inches='tight')

inset_plot_3band(1)
