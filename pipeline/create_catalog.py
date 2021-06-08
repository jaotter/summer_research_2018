from gaussfit_catalog import gaussfit_catalog
import regions
import radio_beam

from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, join, Column
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import scipy.special as special

from scipy.stats import median_abs_deviation
import numpy as np
import os

#this file will create a catalog from n images for each band of data.
#columns (where n is the image name)
#srcID, fwhm_maj_n, fwhm_maj_err_n, fwhm_min_n, fwhm_min_err_n, fwhm_maj_deconv_n, fwhm_maj_deconv_err_n, fwhm_min_deconv_n, fwhm_min_deconv_err_n, aspect_ratio_deconv, aspect_ratio_deconv_err, pa_n, pa_err_n, ap_flux_n, ap_flux_err_n, RA_n, RA_err_n, DEC_n, DEC_err_n


def mask(reg, cutout):#masks everything except the region                                                                                                                                     
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')


def single_img_catalog(B3_img, B3_name, B6_img, B6_name, B7_img, B7_name, cat_name, nonconv_B6_img=None, nonconv_B7_img=None):
    #creates catalog from one image in each band
    #B3_names, B6_names, B7_names only used for gaussian diag directory names
    
    #ref_data_name = '/home/jotter/nrao/summer_research_2018/tables/ref_catalog_added.fits'
    ref_data_name = '/home/jotter/nrao/summer_research_2018/tables/ref_catalog_feb21_updt.fits'
    ref_data = Table.read(ref_data_name)
    ref_arrs = [ref_data['B3_detect'], ref_data['B6_detect'], ref_data['B7_detect']]
  
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

        if name == 'B6' and nonconv_B6_img is not None:
            fl = fits.open(nonconv_B6_img)
            header = fl[0].header
            nonconv_beam = radio_beam.Beam.from_fits_header(header)
            ppbeam = (nonconv_beam.sr/(pixel_scale**2)).decompose().value
        if name == 'B7' and nonconv_B7_img is not None:
            fl = fits.open(nonconv_B7_img)
            header = fl[0].header
            nonconv_beam = radio_beam.Beam.from_fits_header(header)
            ppbeam = (nonconv_beam.sr/(pixel_scale**2)).decompose().value
            
        #now get ready to fit gaussians
        #start by setting up save directory for images
        gauss_save_dir = '/home/jotter/nrao/gauss_diags_feb21/'+img_name+'/'
        if not os.path.exists(gauss_save_dir):
            os.makedirs(gauss_save_dir)
        #now make region list

        rad = Angle(1, 'arcsecond') #radius used in region list
        regs = []

        src_inds = np.where(ref_arrs[b] == True)[0]
        print(len(src_inds))

        for ind in src_inds:
            reg = regions.CircleSkyRegion(center=SkyCoord(ref_data['RA_B3'][ind]*u.degree, ref_data['DEC_B3'][ind]*u.degree), radius=rad, meta={'text':str(ref_data['D_ID'][ind])})
            reg_pix = reg.to_pixel(img_wcs)
            if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                    if np.isnan(img_data[int(reg_pix.center.x), int(reg_pix.center.y)]) == False:
                        regs.append(reg)

        cat_r = Angle(0.5, 'arcsecond')/2 #radius for gaussian fitting
        print('ok')
        gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=gauss_save_dir, max_offset_in_beams = 1, max_radius_in_beams = 5)
        #table does not have all columns yet, add others later
        img_table = Table(names=('D_ID', 'fwhm_maj_'+name, 'fwhm_maj_err_'+name, 'fwhm_min_'+name, 'fwhm_min_err_'+name, 'pa_'+name, 'pa_err_'+name, 'gauss_amp_'+name, 'gauss_amp_err_'+name,'RA_'+name,'RA_err_'+name, 'DEC_'+name, 'DEC_err_'+name, 'success_'+name), dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'S8'))
        for key in gauss_cat:
            img_table.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], gauss_cat[key]['center_x'], gauss_cat[key]['e_center_x'], gauss_cat[key]['center_y'], gauss_cat[key]['e_center_y'], gauss_cat[key]['success']))
        #now measure deconvovled sizes and aperture flux measurements for each source 
        ap_flux_arr = []
        ap_flux_err_arr = []
        fwhm_maj_deconv_arr = []
        fwhm_maj_deconv_err_arr = []
        fwhm_min_deconv_arr = []
        fwhm_min_deconv_err_arr = []
        pa_deconv_arr = []
        pa_deconv_err_arr = []
        snr_arr = []
        rms_arr = []
 
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

                ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm.value*2, pix_minor_fwhm.value*2, angle=pos_ang)
                ap_mask = ellipse_reg.to_mask()
                cutout_mask = ap_mask.cutout(img_data)
                
                aperture_flux = np.nansum(cutout_mask[ap_mask.data==1])/ppbeam
                npix = len(cutout_mask[ap_mask.data==1])

                #now make annulus for measuring background and error
                annulus_width = 15 #pixels
                annulus_radius = img_table['fwhm_maj_'+name][row]*u.arcsecond#+0.05*u.arcsecond
                annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

                #cutout image
                cutout = Cutout2D(img_data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
                cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

                #define aperture regions for SNR
                innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value)
                outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value+annulus_width)

                #Make masks from aperture regions
                annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

                # Calculate the SNR and aperture flux sums
                pixels_in_annulus = cutout.data[annulus_mask.astype('bool')]
                bg_rms = median_abs_deviation(pixels_in_annulus)
                print(img_table['D_ID'][row])
                print('BG RMS: %f' % (bg_rms))
                ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
                bg_median = np.nanmedian(pixels_in_annulus)

                pix_bg = bg_median*npix/ppbeam

                ap_flux_bgcorrect = aperture_flux - pix_bg
                ap_flux_correct = ap_flux_bgcorrect + ap_flux_bgcorrect*(1 - special.erf(2*np.sqrt(np.log(2)))) #flux correction for summing within 2*fwhm
                
                ap_flux_err_arr.append(ap_bg_rms)
                ap_flux_arr.append(ap_flux_correct)
                snr_arr.append(img_table['gauss_amp_'+name][row]/bg_rms)
                rms_arr.append(bg_rms)
                
        cols = ['ap_flux_'+name, 'ap_flux_err_'+name, 'fwhm_maj_deconv_'+name, 'fwhm_maj_deconv_err_'+name, 'fwhm_min_deconv_'+name, 'fwhm_min_deconv_err_'+name, 'pa_deconv_'+name, 'pa_deconv_err_'+name, 'SNR_'+name, 'RMS_'+name]
        arrs = [ap_flux_arr, ap_flux_err_arr, fwhm_maj_deconv_arr, fwhm_maj_deconv_err_arr, fwhm_min_deconv_arr, fwhm_min_deconv_err_arr, pa_deconv_arr, pa_deconv_err_arr, snr_arr, rms_arr]
        for c in range(len(cols)):
            img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
        img_table.add_column(Column(np.array(fwhm_maj_deconv_arr)/np.array(fwhm_min_deconv_arr)), name='ar_deconv_'+name)
        band_tables.append(img_table) #list of tables for each image
            
    
    B3B6 = join(band_tables[0], band_tables[1], keys='D_ID', join_type='outer')
    all_bands = join(B3B6, band_tables[2], keys='D_ID', join_type='outer')

    all_bands.write('/home/jotter/nrao/summer_research_2018/tables/'+cat_name+'.fits',  overwrite=True)


def huge_B3_catalog(B3_img, cat_name):
    #creates catalog from one image in each band
    #B3_names, B6_names, B7_names only used for gaussian diag directory names
    
    #ref_data_name = '/home/jotter/nrao/summer_research_2018/tables/ref_catalog_added.fits'
    ref_data_name = '/lustre/cv/observers/cv-12578/orion_disks/tables/ref_catalog_may21.fits'
    ref_data = Table.read(ref_data_name)

    band_tables = []

    name = 'B3'
    fl = fits.open(B3_img)
    header = fl[0].header
    img_data = fl[0].data.squeeze()
    img_wcs = WCS(header).celestial

    beam = radio_beam.Beam.from_fits_header(header)
    pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

    #now get ready to fit gaussians
    #start by setting up save directory for images
    gauss_save_dir = '/lustre/cv/observers/cv-12578/orion_disks/gauss_diags_may21/B3_huge/'
    if not os.path.exists(gauss_save_dir):
        os.makedirs(gauss_save_dir)
    #now make region list

    rad = Angle(0.5, 'arcsecond') #radius used in region list
    regs = []

    for ind in range(len(ref_data)):
        reg = regions.CircleSkyRegion(center=SkyCoord(ref_data['RA_B3'][ind]*u.degree, ref_data['DEC_B3'][ind]*u.degree, frame='icrs'), radius=rad, meta={'text':str(ref_data['Seq_B3'][ind])})
        reg_pix = reg.to_pixel(img_wcs)
        if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
            if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                if np.isnan(img_data[int(reg_pix.center.x), int(reg_pix.center.y)]) == False:
                    regs.append(reg)

    cat_r = Angle(0.5, 'arcsecond') #radius for gaussian fitting
    print(len(regs))
    gauss_cat = gaussfit_catalog(B3_img, regs, cat_r, savepath=gauss_save_dir, max_offset_in_beams = 1, max_radius_in_beams = 5)
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
    snr_arr = []
    rms_arr = []

    print(len(img_table))
    for row in range(len(img_table)): #now loop through sources in reference data and make measurements
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

        ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm.value*2, pix_minor_fwhm.value*2, angle=pos_ang)
        ap_mask = ellipse_reg.to_mask()
        cutout_mask = ap_mask.cutout(img_data)
        
        aperture_flux = np.nansum(cutout_mask[ap_mask.data==1])/ppbeam
        npix = len(cutout_mask[ap_mask.data==1])

        #now make annulus for measuring background and error
        annulus_width = 15 #pixels
        annulus_radius = img_table['fwhm_maj_'+name][row]*u.arcsecond#+0.05*u.arcsecond
        annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()
        
        #cutout image
        cutout = Cutout2D(img_data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
        cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

        #define aperture regions for SNR
        innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value)
        outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value+annulus_width)
        
        #Make masks from aperture regions
        annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

        # Calculate the SNR and aperture flux sums
        pixels_in_annulus = cutout.data[annulus_mask.astype('bool')]
        bg_rms = median_abs_deviation(pixels_in_annulus)
        #print(img_table['D_ID'][row])
        #print('BG RMS: %f' % (bg_rms))
        ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
        bg_median = np.nanmedian(pixels_in_annulus)
        
        pix_bg = bg_median*npix/ppbeam
            
        ap_flux_bgcorrect = aperture_flux - pix_bg
        ap_flux_correct = ap_flux_bgcorrect + ap_flux_bgcorrect*(1 - special.erf(2*np.sqrt(np.log(2)))) #flux correction for summing within 2*fwhm

        ap_flux_err_arr.append(ap_bg_rms)
        ap_flux_arr.append(ap_flux_correct)
        snr_arr.append(img_table['gauss_amp_'+name][row]/bg_rms)
        rms_arr.append(bg_rms)

    cols = ['ap_flux_'+name, 'ap_flux_err_'+name, 'fwhm_maj_deconv_'+name, 'fwhm_maj_deconv_err_'+name, 'fwhm_min_deconv_'+name, 'fwhm_min_deconv_err_'+name, 'pa_deconv_'+name, 'pa_deconv_err_'+name, 'SNR_'+name, 'RMS_'+name]
    arrs = [ap_flux_arr, ap_flux_err_arr, fwhm_maj_deconv_arr, fwhm_maj_deconv_err_arr, fwhm_min_deconv_arr, fwhm_min_deconv_err_arr, pa_deconv_arr, pa_deconv_err_arr, snr_arr, rms_arr]
    for c in range(len(cols)):
        print(len(arrs[c]))
        img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
    img_table.add_column(Column(np.array(fwhm_maj_deconv_arr)/np.array(fwhm_min_deconv_arr)), name='ar_deconv_'+name)
    img_table.write('/lustre/cv/observers/cv-12578/orion_disks/tables/'+cat_name+'.fits',  overwrite=True)



def b6b7_catalog(B6_img, B6_name, B7_img, B7_name, cat_name, nonconv_B6_img=None, nonconv_B7_img=None):
    #creates catalog from one image in each band
    #B3_names, B6_names, B7_names only used for gaussian diag directory names
    
    ref_data_name = '/home/jotter/nrao/summer_research_2018/tables/ref_catalog_may21.fits'
    ref_data = Table.read(ref_data_name)
    ref_arrs = [ref_data['B6_detect'], ref_data['B7_detect']]
  
    band_imgs = [B6_img, B7_img]
    band_names = ['B6', 'B7']
    band_img_names = [B6_name, B7_name]
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

        if name == 'B6' and nonconv_B6_img is not None:
            fl = fits.open(nonconv_B6_img)
            header = fl[0].header
            nonconv_beam = radio_beam.Beam.from_fits_header(header)
            ppbeam = (nonconv_beam.sr/(pixel_scale**2)).decompose().value
        if name == 'B7' and nonconv_B7_img is not None:
            fl = fits.open(nonconv_B7_img)
            header = fl[0].header
            nonconv_beam = radio_beam.Beam.from_fits_header(header)
            ppbeam = (nonconv_beam.sr/(pixel_scale**2)).decompose().value
            
        #now get ready to fit gaussians
        #start by setting up save directory for images
        gauss_save_dir = '/home/jotter/nrao/gauss_diags_may21/'+img_name+'/'
        if not os.path.exists(gauss_save_dir):
            os.makedirs(gauss_save_dir)
        #now make region list

        rad = Angle(1, 'arcsecond') #radius used in region list
        regs = []

        src_inds = np.where(ref_arrs[b] == True)[0]
        print(len(src_inds))

        for ind in src_inds:
            reg = regions.CircleSkyRegion(center=SkyCoord(ref_data['RA_B3'][ind]*u.degree, ref_data['DEC_B3'][ind]*u.degree), radius=rad, meta={'text':str(ref_data['B3_Seq'][ind])})
            reg_pix = reg.to_pixel(img_wcs)
            if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                    if np.isnan(img_data[int(reg_pix.center.x), int(reg_pix.center.y)]) == False:
                        regs.append(reg)

        cat_r = Angle(0.5, 'arcsecond')/2 #radius for gaussian fitting
        print('ok')
        gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=gauss_save_dir, max_offset_in_beams = 1, max_radius_in_beams = 5)
        #table does not have all columns yet, add others later
        img_table = Table(names=('Seq_B3', 'fwhm_maj_'+name, 'fwhm_maj_err_'+name, 'fwhm_min_'+name, 'fwhm_min_err_'+name, 'pa_'+name, 'pa_err_'+name, 'gauss_amp_'+name, 'gauss_amp_err_'+name,'RA_'+name,'RA_err_'+name, 'DEC_'+name, 'DEC_err_'+name, ), dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
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
        snr_arr = []
        rms_arr = []
 
        for row in range(len(img_table)): #now loop through sources in reference data and make measurements
            ref_ind = np.where(ref_data['B3_Seq'] == img_table['Seq_B3'][row])[0]
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

                ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm.value*2, pix_minor_fwhm.value*2, angle=pos_ang)
                ap_mask = ellipse_reg.to_mask()
                cutout_mask = ap_mask.cutout(img_data)
                
                aperture_flux = np.nansum(cutout_mask[ap_mask.data==1])/ppbeam
                npix = len(cutout_mask[ap_mask.data==1])

                #now make annulus for measuring background and error
                annulus_width = 15 #pixels
                annulus_radius = img_table['fwhm_maj_'+name][row]*u.arcsecond#+0.05*u.arcsecond
                annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

                #cutout image
                cutout = Cutout2D(img_data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
                cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

                #define aperture regions for SNR
                innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value)
                outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value+annulus_width)

                #Make masks from aperture regions
                annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

                # Calculate the SNR and aperture flux sums
                pixels_in_annulus = cutout.data[annulus_mask.astype('bool')]
                bg_rms = median_abs_deviation(pixels_in_annulus)
                print(img_table['Seq_B3'][row])
                print('BG RMS: %f' % (bg_rms))
                ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
                bg_median = np.nanmedian(pixels_in_annulus)

                pix_bg = bg_median*npix/ppbeam

                ap_flux_bgcorrect = aperture_flux - pix_bg
                ap_flux_correct = ap_flux_bgcorrect + ap_flux_bgcorrect*(1 - special.erf(2*np.sqrt(np.log(2)))) #flux correction for summing within 2*fwhm
                
                ap_flux_err_arr.append(ap_bg_rms)
                ap_flux_arr.append(ap_flux_correct)
                snr_arr.append(img_table['gauss_amp_'+name][row]/bg_rms)
                rms_arr.append(bg_rms)
                
        cols = ['ap_flux_'+name, 'ap_flux_err_'+name, 'fwhm_maj_deconv_'+name, 'fwhm_maj_deconv_err_'+name, 'fwhm_min_deconv_'+name, 'fwhm_min_deconv_err_'+name, 'pa_deconv_'+name, 'pa_deconv_err_'+name, 'SNR_'+name, 'RMS_'+name]
        arrs = [ap_flux_arr, ap_flux_err_arr, fwhm_maj_deconv_arr, fwhm_maj_deconv_err_arr, fwhm_min_deconv_arr, fwhm_min_deconv_err_arr, pa_deconv_arr, pa_deconv_err_arr, snr_arr, rms_arr]
        for c in range(len(cols)):
            img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
        img_table.add_column(Column(np.array(fwhm_maj_deconv_arr)/np.array(fwhm_min_deconv_arr)), name='ar_deconv_'+name)
        band_tables.append(img_table) #list of tables for each image
            
    
    B6B7 = join(band_tables[0], band_tables[1], keys='Seq_B3', join_type='outer')
    B6B7.write('/home/jotter/nrao/summer_research_2018/tables/'+cat_name+'.fits',  overwrite=True)
