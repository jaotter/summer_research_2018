from gaussfit_catalog import gaussfit_catalog
import regions
import radio_beam

from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, join, Column
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

import numpy as np
import os
from dendrogram_catalog import mask, rms

#this file will create a catalog from n images for each band of data.
#columns (where n is the image name)
#srcID, fwhm_maj_n, fwhm_maj_err_n, fwhm_min_n, fwhm_min_err_n, fwhm_maj_deconv_n, fwhm_maj_deconv_err_n, fwhm_min_deconv_n, fwhm_min_deconv_err_n, aspect_ratio_deconv, aspect_ratio_deconv_err, pa_n, pa_err_n, ap_flux_n, ap_flux_err_n, RA_n, RA_err_n, DEC_n, DEC_err_n


def single_img_catalog(B3_img, B3_name, B6_img, B6_name, B7_img, B7_name, cat_name):
    #creates catalog from one image in each band
    #B3_names, B6_names, B7_names only used for gaussian diag directory names
    
    ref_data_name = '/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits' #master_500klplus_B3_ref.fits'
    ref_data = Table.read(ref_data_name)
    ref_arrs = [ref_data['_idx_B3'], ref_data['_idx_B6'], ref_data['_idx_B7_hr']]
  
    band_imgs = [B3_img, B6_img, B7_img]
    band_names = ['B3', 'B6', 'B7']
    band_img_names = [B3_name, B6_name, B7_name]
    band_tables = []
    for b in range(len(band_imgs)): #first loop through different bands
        name = band_names[b]
        img_name = band_img_names[b]
        img = band_imgs[b]
        fl = fits.open(img)
        header = fl[0].header
        img_data = fl[0].data.squeeze()
        img_wcs = WCS(header).celestial

        beam = radio_beam.Beam.from_fits_header(header)
        pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
        if name == 'B6':
            ppbeam = 127
        if name == 'B7':
            ppbeam = 51
        #now get ready to fit gaussians
        #start by setting up save directory for images
        gauss_save_dir = '/lustre/aoc/students/jotter/gauss_diags/create_cat/'+img_name+'/'
        if not os.path.exists(gauss_save_dir):
            os.makedirs(gauss_save_dir)
        #now make region list

        rad = Angle(1, 'arcsecond') #radius used in region list
        regs = []

        src_inds = np.where(np.isnan(ref_arrs[b]) == False)[0]
        print(len(src_inds))

        band_ind_nm = name
        if name == 'B7':
            band_ind_nm = 'B7_hr'
        for ind in src_inds:
            reg = regions.CircleSkyRegion(center=SkyCoord(ref_data['gauss_x_'+band_ind_nm][ind]*u.degree, ref_data['gauss_y_'+band_ind_nm][ind]*u.degree), radius=rad, meta={'text':str(ref_data['D_ID'][ind])})
            reg_pix = reg.to_pixel(img_wcs)
            if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                    if np.isnan(img_data[int(reg_pix.center.x), int(reg_pix.center.y)]) == False:
                        regs.append(reg)

        cat_r = Angle(0.5, 'arcsecond') #radius for gaussian fitting
        gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=gauss_save_dir, max_offset_in_beams = 1, max_radius_in_beams = 5)
        #table does not have all columns yet, add others later
        img_table = Table(names=('D_ID', 'fwhm_maj_'+name, 'fwhm_maj_err_'+name, 'fwhm_min_'+name, 'fwhm_min_err_'+name, 'pa_'+name, 'pa_err_'+name, 'gauss_amp_'+name, 'gauss_amp_err_'+name,'RA_'+name,'RA_err_'+name, 'DEC_'+name, 'DEC_err_'+name), dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
        for key in gauss_cat:
            img_table.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], gauss_cat[key]['center_x'], gauss_cat[key]['e_center_x'], gauss_cat[key]['center_y'], gauss_cat[key]['e_center_y']))
        #now measure deconvovled sizes and aperture flux measurements for each source 
        ap_flux_arr = []
        ap_flux_err_arr = []
        fwhm_maj_deconv_arr = []
        fwhm_maj_deconv_err_arr = []
        fwhm_min_deconv_arr = []
        fwhm_min_deconv_err_arr = []
        pa_deconv_arr = []
        pa_deconv_err_arr = []
 
        for row in range(len(img_table)): #now loop through sources in reference data and make measurements
            ref_ind = np.where(ref_data['D_ID'] == img_table['D_ID'][row])[0]
            if len(ref_ind > 0):
                #now measuring deconvolved sizes
                measured_source_size = radio_beam.Beam(major=img_table['fwhm_maj_'+name][row]*u.arcsec, minor=img_table['fwhm_min_'+name][row]*u.arcsec, pa=(img_table['pa_'+name][row]-90)*u.degree)
                try:
                    deconv_size = measured_source_size.deconvolve(beam)
                    fwhm_maj_deconv_arr.append(deconv_size.major.value)
                    fwhm_min_deconv_arr.append(deconv_size.minor.value)
                    fwhm_maj_deconv_err_arr.append(img_table['fwhm_maj_err_'+name][row]) #same error as non deconvolved
                    fwhm_min_deconv_err_arr.append(img_table['fwhm_min_err_'+name][row])
                    pa_deconv_arr.append(deconv_size.pa.to(u.deg).value)
                    pa_deconv_err_arr.append(img_table['pa_err_'+name][row])
                except ValueError:
                    fwhm_maj_deconv_arr.append(np.nan)
                    fwhm_min_deconv_arr.append(np.nan)
                    fwhm_maj_deconv_err_arr.append(np.nan)
                    fwhm_min_deconv_err_arr.append(np.nan)
                    pa_deconv_arr.append(np.nan)
                    pa_deconv_err_arr.append(np.nan)
 

                pix_major_fwhm = ((img_table['fwhm_maj_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                pix_minor_fwhm = ((img_table['fwhm_min_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                center_coord = SkyCoord(img_table['RA_'+name][row], img_table['DEC_'+name][row], frame='icrs', unit=(u.deg, u.deg))
                center_coord_pix = center_coord.to_pixel(img_wcs)
                center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])
                pos_ang = (img_table['pa_'+name][row]-90)*u.deg #must subtract 90 to be consistent
                ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm*2, pix_minor_fwhm*2, angle=pos_ang)
                ap_mask = ellipse_reg.to_mask()
                cutout_mask = ap_mask.cutout(img_data)
                
                aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
                npix = len(cutout_mask[ap_mask.data==1])

                #now make annulus for measuring background and error
                annulus_width = 15 #pixels
                annulus_radius = 0.1*u.arcsecond
                annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

                #cutout image
                cutout = Cutout2D(img_data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
                cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

                #define aperture regions for SNR
                innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix)
                outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix+annulus_width)

                #Make masks from aperture regions
                annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

                # Calculate the SNR and aperture flux sums
                pixels_in_annulus = cutout.data[annulus_mask.astype('bool')]
                bg_rms = rms(pixels_in_annulus)
                ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
                bg_median = np.median(pixels_in_annulus)

                pix_bg = bg_median*npix/ppbeam
                
                ap_flux_err_arr.append(ap_bg_rms)
                ap_flux_arr.append(aperture_flux - pix_bg)
                
        cols = ['ap_flux_'+name, 'ap_flux_err_'+name, 'fwhm_maj_deconv_'+name, 'fwhm_maj_deconv_err_'+name, 'fwhm_min_deconv_'+name, 'fwhm_min_deconv_err_'+name, 'pa_deconv_'+name, 'pa_deconv_err_'+name]
        arrs = [ap_flux_arr, ap_flux_err_arr, fwhm_maj_deconv_arr, fwhm_maj_deconv_err_arr, fwhm_min_deconv_arr, fwhm_min_deconv_err_arr, pa_deconv_arr, pa_deconv_err_arr]
        for c in range(len(cols)):
            img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
        img_table.add_column(Column(np.array(fwhm_maj_deconv_arr)/np.array(fwhm_min_deconv_arr)), name='ar_deconv_'+name)
        band_tables.append(img_table) #list of tables for each image
            
    
    B3B6 = join(band_tables[0], band_tables[1], keys='D_ID', join_type='outer')
    all_bands = join(B3B6, band_tables[2], keys='D_ID', join_type='outer')

    all_bands.write('/users/jotter/summer_research_2018/tables/'+cat_name+'.fits',  overwrite=True)





















def imgs_catalog(B3_imgs, B3_names, B6_imgs, B6_names, B7_imgs, B7_names, short_names, cat_name):

    #B3_names, B6_names, B7_names only used for gaussian diag directory names
    #short names actually used in table
    #creates catalog from multiple images in each band
    
    ref_data_name = '/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits' #master_500klplus_B3_ref.fits'
    ref_data = Table.read(ref_data_name)
    ref_arrs = [ref_data['_idx_B3'], ref_data['_idx_B6'], ref_data['_idx_B7_hr']]
    
    band_imgs = [B3_imgs, B6_imgs, B7_imgs]
    band_names = ['B3', 'B6', 'B7']
    band_img_names = [B3_names, B6_names, B7_names]
    band_tables = []
    for b in range(len(band_imgs)): #first loop through different bands
        band = band_imgs[b]
        band_nm = band_names[b]
        img_names = band_img_names[b]
        img_tables = []
        for img_n in range(len(band)): #now loop through images for each band
            file_name = img_names[img_n]
            name = band_nm+'_'+short_names[img_n]
            img = band[img_n]
            fl = fits.open(img)
            header = fl[0].header
            img_data = fl[0].data.squeeze()
            img_wcs = WCS(header).celestial

            beam = radio_beam.Beam.from_fits_header(header)
            pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
            ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

            #now get ready to fit gaussians
            #start by setting up save directory for images
            gauss_save_dir = '/lustre/aoc/students/jotter/gauss_diags/create_cat/'+file_name+'/'
            if not os.path.exists(gauss_save_dir):
                os.makedirs(gauss_save_dir)
            #now make region list

            rad = Angle(1, 'arcsecond') #radius used in region list
            regs = []

            src_inds = np.where(np.isnan(ref_arrs[b]) == False)[0]
            print(len(src_inds))

            band_ind_nm = band_nm
            if band_nm == 'B7':
                band_ind_nm = 'B7_hr'
            for ind in src_inds:
                reg = regions.CircleSkyRegion(center=SkyCoord(ref_data['gauss_x_'+band_ind_nm][ind]*u.degree, ref_data['gauss_y_'+band_ind_nm][ind]*u.degree), radius=rad, meta={'text':str(ref_data['D_ID'][ind])})
                reg_pix = reg.to_pixel(img_wcs)
                if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                    if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                        if np.isnan(img_data[int(reg_pix.center.x), int(reg_pix.center.y)]) == False:
                            regs.append(reg)

            cat_r = Angle(0.5, 'arcsecond') #radius for gaussian fitting
            gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=gauss_save_dir)#, max_radius_in_beams = 15)
            #table does not have all columns yet, add others later
            img_table = Table(names=('D_ID', 'fwhm_maj_'+name, 'fwhm_maj_err_'+name, 'fwhm_min_'+name, 'fwhm_min_err_'+name, 'pa_'+name, 'pa_err_'+name, 'gauss_amp_'+name, 'gauss_amp_err_'+name, 'RA_'+name,'RA_err_'+name, 'DEC_'+name, 'DEC_err_'+name), dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
            for key in gauss_cat:
                img_table.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], gauss_cat[key]['center_x'], gauss_cat[key]['e_center_x'], gauss_cat[key]['center_y'], gauss_cat[key]['e_center_y']))
            #now measure deconvovled sizes and aperture flux measurements for each source 
            ap_flux_arr = []
            ap_flux_err_arr = []
            fwhm_maj_deconv_arr = []
            fwhm_maj_deconv_err_arr = []
            fwhm_min_deconv_arr = []
            fwhm_min_deconv_err_arr = []
            ar_deconv_arr = []
            ar_deconv_err_arr = []
            
            for row in range(len(img_table)): #now loop through sources in reference data and make measurements
                ref_ind = np.where(ref_data['D_ID'] == img_table['D_ID'][row])[0]
                if len(ref_ind > 0):
                    pix_major_fwhm = ((img_table['fwhm_maj_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                    pix_minor_fwhm = ((img_table['fwhm_min_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                    center_coord = SkyCoord(img_table['RA_'+name][row], img_table['DEC_'+name][row], frame='icrs', unit=(u.deg, u.deg))
                    center_coord_pix = center_coord.to_pixel(img_wcs)
                    center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])

                    pos_ang = img_table['pa_'+name][row]*u.deg

                    ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm*2, pix_minor_fwhm*2, angle=pos_ang)
                    size = pix_major_fwhm*2.1
                    ap_mask = ellipse_reg.to_mask()
                    cutout_mask = ap_mask.cutout(img_data)

                    aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
                    npix = len(cutout_mask[ap_mask.data==1])
                    
                    #now make annulus for measuring background and error
                    annulus_width = 15 #pixels
                    annulus_radius = 0.1*u.arcsecond
                    annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()
                    
                    #cutout image
                    cutout = Cutout2D(img_data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
                    cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

                    #define aperture regions for SNR
                    innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix)
                    outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix+annulus_width)

                    #Make masks from aperture regions
                    annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

                    # Calculate the SNR and aperture flux sums
                    pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
                    bg_rms = rms(pixels_in_annulus)
                    ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
                    bg_median = np.median(pixels_in_annulus)

                    pix_bg = bg_median*npix/ppbeam

                    ap_flux_err_arr.append(ap_bg_rms)
                    ap_flux_arr.append(aperture_flux - pix_bg)

                    #now measuring deconvolved sizes
                    measured_source_size = radio_beam.Beam(major=img_table['fwhm_maj_'+name][row]*u.arcsec, minor=img_table['fwhm_min_'+name][row]*u.arcsec, pa=(img_table['pa_'+name][row]-90)*u.degree)
                    try:
                        deconv_size = measured_source_size.deconvolve(beam)
                        fwhm_maj_deconv_arr.append(deconv_size.major.value)
                        fwhm_min_deconv_arr.append(deconv_size.minor.value)
                        fwhm_maj_deconv_err_arr.append(img_table['fwhm_maj_err_'+name][row]) #same error as non deconvolved
                        fwhm_min_deconv_err_arr.append(img_table['fwhm_min_err_'+name][row])
                        #aspect_ratio = deconv_size.major.value/deconv_size.minor.value
                        #ar_deconv_arr.append(aspect_ratio)
                        #ar_err = np.sqrt((img_table['fwhm_maj_err_'+name]/deconv_size.major.value)**2 + (img_table['fwhm_min_err_'+name]/deconv_si
                    except ValueError:
                        fwhm_maj_deconv_arr.append(np.nan)
                        fwhm_min_deconv_arr.append(np.nan)
                        fwhm_maj_deconv_err_arr.append(np.nan)
                        fwhm_min_deconv_err_arr.append(np.nan)
                        #pa_deconv_arr.append(np.nan)
                        #pa_deconv_err_arr.append(np.nan)
                        
            cols = ['ap_flux_'+name, 'ap_flux_err_'+name, 'fwhm_maj_deconv_'+name, 'fwhm_maj_deconv_err_'+name, 'fwhm_min_deconv_'+name, 'fwhm_min_deconv_err_'+name]#, 'ar_deconv_'+name, 'ar_deconv_err_'+name]
            arrs = [ap_flux_arr, ap_flux_err_arr, fwhm_maj_deconv_arr, fwhm_maj_deconv_err_arr, fwhm_min_deconv_arr, fwhm_min_deconv_err_arr]#, ar_deconv_arr, ar_deconv_err_arr]
            for c in range(len(cols)):
                img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
            img_table.add_column(Column(np.array(fwhm_maj_deconv_arr)/np.array(fwhm_min_deconv_arr)), name='ar_deconv_'+name)
            img_tables.append(img_table) #list of tables for each image
            
            
        band_table = join(img_tables[0], img_tables[1], keys='D_ID', join_type='outer')#table combining all images for a band of data
        band_tables.append(band_table)
    B3B6 = join(band_tables[0], band_tables[1], keys='D_ID', join_type='outer')
    all_bands = join(B3B6, band_tables[2], keys='D_ID', join_type='outer')

    all_bands.write('/lustre/aoc/students/jotter/dendro_catalogs/'+cat_name+'.fits',  overwrite=True)
    #lastly, match with other data
