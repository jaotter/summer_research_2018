from multi_gauss_fit import bg_gaussfit
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

#this file fits a source and background with two gaussians
#columns (where n is the image name)
#srcID, fwhm_maj_n, fwhm_maj_err_n, fwhm_min_n, fwhm_min_err_n, fwhm_maj_deconv_n, fwhm_maj_deconv_err_n, fwhm_min_deconv_n, fwhm_min_deconv_err_n, aspect_ratio_deconv, aspect_ratio_deconv_err, pa_n, pa_err_n, ap_flux_n, ap_flux_err_n, RA_n, RA_err_n, DEC_n, DEC_err_n


def fit_source(srcID, img, img_name, band, bg_stddev_x, bg_stddev_y):
    #this function fits a given source, and the background
    #srcID : int
        #name of source to fit in catalogs
    #img : fits file
        #fits file with source to fit
    #img_name : str
        #name of image for the directory where the fit plots will go
    #band : str
        #band of image to fit ('B3', 'B6', or 'B7')
    #bg_stddev_x : float
        #eyeballed estimate of stddev of the background source in pixels
    #bg_stddev_y : float
        #same as above in y direction
    
    ref_data_name = '/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_flags.fits'
    ref_data = Table.read(ref_data_name)
    
    fl = fits.open(img)
    header = fl[0].header
    img_data = fl[0].data.squeeze()
    img_wcs = WCS(header).celestial

    beam = radio_beam.Beam.from_fits_header(header)
    pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

    #now get ready to fit gaussians
    #start by setting up save directory for images
    gauss_save_dir = '/lustre/aoc/students/jotter/gauss_diags/create_cat/'+img_name+'/'
    if not os.path.exists(gauss_save_dir):
        os.makedirs(gauss_save_dir)
    #now make region
    rad = Angle(1, 'arcsecond') #radius used in region list
    src_ind = np.where(ref_data['D_ID']==srcID)[0]
    
    reg = regions.CircleSkyRegion(center=SkyCoord(ref_data['RA_'+band][src_ind], ref_data['DEC_'+band][src_ind], unit='deg'), radius=rad, meta={'text':str(ref_data['D_ID'][src_ind])+'xstddev_'+str(bg_stddev_x)+'_ystddev_'+str(bg_stddev_y)})

    cat_r = Angle(0.5, 'arcsecond') #radius for gaussian fitting
    gauss_cat = bg_gaussfit(img, reg, cat_r, bg_stddev_x=bg_stddev_x, bg_stddev_y=bg_stddev_y, savepath=gauss_save_dir, max_offset_in_beams = 1, max_radius_in_beams = 5)

    img_table = Table(names=('D_ID', 'fwhm_maj_'+name, 'fwhm_maj_err_'+name, 'fwhm_min_'+name, 'fwhm_min_err_'+name, 'pa_'+name, 'pa_err_'+name, 'fwhm_maj_deconv_'+name, 'fwhm_min_deconv_'+name, 'deconv_pa_'+name, 'RA_'+name,'RA_err_'+name, 'DEC_'+name, 'DEC_err_'+name), dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
    for key in gauss_cat:
        img_table.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['deconv_fwhm_major'], gauss_cat[key]['deconv_fwhm_minor'], gauss_cat[key]['deconv_pa'], gauss_cat[key]['center_x'], gauss_cat[key]['e_center_x'], gauss_cat[key]['center_y'], gauss_cat[key]['e_center_y']))

    #now measure deconvovled sizes and aperture flux measurements for each source 
    ap_flux_arr = []
    ap_flux_err_arr = []

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

    cols = ['ap_flux_'+band, 'ap_flux_err_'+band]
    arrs = [ap_flux_arr, ap_flux_err_arr]
    for c in range(len(cols)):
        img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
    img_table.add_column(Column(img_table['fwhm_maj_deconv'+band]/img_table['fwhm_min_deconv'+band]), name='ar_deconv_'+band)
    
    return img_table
