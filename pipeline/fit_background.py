from multi_gauss_fit import bg_gaussfit, gaussfit_cutoutim
from gaussfit_catalog import gaussfit_catalog
import regions
import radio_beam

from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, join, Column
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from scipy.stats import median_abs_deviation
import scipy.special as special
import numpy as np
import os

#this file fits a source and background with two gaussians
#columns (where n is the image name)
#srcID, fwhm_maj_n, fwhm_maj_err_n, fwhm_min_n, fwhm_min_err_n, fwhm_maj_deconv_n, fwhm_maj_deconv_err_n, fwhm_min_deconv_n, fwhm_min_deconv_err_n, aspect_ratio_deconv, aspect_ratio_deconv_err, pa_n, pa_err_n, ap_flux_n, ap_flux_err_n, RA_n, RA_err_n, DEC_n, DEC_err_n

def mask(reg, cutout):#masks everything except the region                                                                                                                                     
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')


def fit_source(srcID, img, img_name, band, fit_bg=False, bg_stddev_x=30, bg_stddev_y=30, bg_mean_x=0, bg_mean_y=0, zoom=1, max_offset_in_beams=1, max_radius_in_beams=5, nonconv_img=None, mask_size=1.5):
    #this function fits a given source, and the background
    #srcID : int
        #name of source to fit in catalogs
    #img : fits file
        #fits file with source to fit
    #img_name : str
        #name of image for the directory where the fit plots will go
    #band : str
        #band of image to fit ('B3', 'B6', or 'B7')
    #fit_bg : bool
    #if False, do not fit background gaussian
    #bg_stddev_x : float
        #eyeballed estimate of stddev of the background source in pixels
    #bg_stddev_y : float
        #same as above in y direction
    #bg_mean_x/y : float
        #pixels away from center (origin) in x/y direction for background gaussian mean guess
    #zoom : float
        #amount of zoom, values greater than 1 are zoom ins

    
    #ref_data_name = '/home/jotter/nrao/summer_research_2018/tables/dendro_ref_catalog_edited.fits'
    #ref_data_name = '/home/jotter/nrao/summer_research_2018/tables/ref_catalog_may21.fits'
    ref_data_name = '/lustre/cv/observers/cv-12578/orion_disks/summer_research_2018/tables/ref_catalog_may21.fits'
    ref_data = Table.read(ref_data_name)
    
    fl = fits.open(img)
    header = fl[0].header
    img_data = fl[0].data.squeeze()
    img_wcs = WCS(header).celestial

    beam = radio_beam.Beam.from_fits_header(header)
    pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

    if nonconv_img is not None:
        flnonconv = fits.open(nonconv_img)
        nonconv_header = flnonconv[0].header
        nonconv_beam = radio_beam.Beam.from_fits_header(nonconv_header)
        ppbeam = (nonconv_beam.sr/(pixel_scale**2)).decompose().value
        
    #now get ready to fit gaussians
    #start by setting up save directory for images
    #gauss_save_dir = '/home/jotter/nrao/gauss_diags_may21/fitbg/'+img_name+'/'
    gauss_save_dir = f'/lustre/cv/observers/cv-12578/orion_disks/gauss_diags_may21/{img_name}/'
    
    print('saving plots to '+gauss_save_dir)
    if not os.path.exists(gauss_save_dir):
        os.makedirs(gauss_save_dir)
    #now make region
    rad = Angle(1, 'arcsecond') #radius used in region list
    #src_ind = np.where(ref_data['D_ID']==srcID)[0]
    src_ind = np.where(ref_data['B3_Seq']==srcID)[0]

    ra = ref_data['RA_B3'][src_ind].data[0]
    dec = ref_data['DEC_B3'][src_ind].data[0]
    center_reg = SkyCoord(ra, dec, unit='deg', frame='icrs')
    reg = regions.CircleSkyRegion(center=center_reg, radius=1*u.arcsecond, meta={'text':str(ref_data['B3_Seq'][src_ind].data[0])+'_xstddev_'+str(bg_stddev_x)+'_ystddev_'+str(bg_stddev_y)})
        
    region_list = []
    #valid_inds = np.where(np.isnan(ref_data[band+'_detect']) == False)[0]
    for ind in range(len(ref_data)):#valid_inds:
        if ref_data['B3_Seq'][ind] == srcID:
            continue
        ra_i = ref_data['RA_B3'][ind]
        dec_i = ref_data['DEC_B3'][ind]
        region_i = regions.CircleSkyRegion(center=SkyCoord(ra_i, dec_i, unit='deg', frame='icrs'), radius=1*u.arcsecond)
        region_list.append(region_i)

    #print(region_list)
    #print(reg)
        
    cat_r = Angle(0.5, 'arcsecond')/zoom #radius for gaussian fitting
    
    if fit_bg == True:
        gauss_cat, fitim_bg = bg_gaussfit(img, reg, region_list, cat_r, bg_stddev_x=bg_stddev_x, bg_stddev_y=bg_stddev_y, bg_mean_x=bg_mean_x, bg_mean_y=bg_mean_y, savepath=gauss_save_dir, max_offset_in_beams = max_offset_in_beams, max_offset_in_beams_bg = 10, max_radius_in_beams = max_radius_in_beams, mask_size=mask_size)

        #print('gauss_cat length ',len(gauss_cat))
        #k = list(gauss_cat.keys())[0]
        #if gauss_cat[k]['success'] == False:
        #    gauss_cat = gaussfit_cutoutim(img, fitim_bg, reg, region_list, cat_r, savepath=gauss_save_dir, max_offset_in_beams = max_offset_in_beams, max_radius_in_beams = max_radius_in_beams)
        #    success = gauss_cat[k]['success']
        #    print(F'ALTERNATIVE FIT SUCCESS: {success}')

    else:
        gauss_cat = gaussfit_catalog(img, [reg], cat_r, savepath=gauss_save_dir, max_offset_in_beams = max_offset_in_beams, max_radius_in_beams = max_radius_in_beams)

    img_table = Table(names=('Seq', 'fwhm_maj_'+band, 'fwhm_maj_err_'+band, 'fwhm_min_'+band, 'fwhm_min_err_'+band, 'pa_'+band, 'pa_err_'+band, 'gauss_amp_'+band, 'gauss_amp_err_'+band, 'RA_'+band,'RA_err_'+band, 'DEC_'+band, 'DEC_err_'+band), dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
    for key in gauss_cat:
        img_table.add_row((srcID, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], gauss_cat[key]['center_x'], gauss_cat[key]['e_center_x'], gauss_cat[key]['center_y'], gauss_cat[key]['e_center_y']))

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
    
    for row in range(len(img_table)): #now loop through sources in reference data and make measurements
        ref_ind = np.where(ref_data['B3_Seq'] == img_table['Seq'][row])[0]
        if True==True:#len(ref_ind > 0):

            measured_source_size = radio_beam.Beam(major=img_table['fwhm_maj_'+band][row]*u.arcsec, minor=img_table['fwhm_min_'+band][row]*u.arcsec, pa=(img_table['pa_'+band][row]-90)*u.deg)
            try:
                deconv_size = measured_source_size.deconvolve(beam)
                fwhm_maj_deconv_arr.append(deconv_size.major.value)
                fwhm_min_deconv_arr.append(deconv_size.minor.value)
                fwhm_maj_deconv_err_arr.append(img_table['fwhm_maj_err_'+band][row])
                fwhm_min_deconv_err_arr.append(img_table['fwhm_min_err_'+band][row])
                pa_deconv_arr.append(deconv_size.pa.value)
                pa_deconv_err_arr.append(img_table['pa_err_'+band][row])

                
            except ValueError:
                fwhm_maj_deconv_arr.append(np.nan)
                fwhm_min_deconv_arr.append(np.nan)
                fwhm_maj_deconv_err_arr.append(np.nan)
                fwhm_min_deconv_err_arr.append(np.nan)
                pa_deconv_arr.append(np.nan)
                pa_deconv_err_arr.append(np.nan)



            pix_major_fwhm = ((img_table['fwhm_maj_'+band][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
            pix_minor_fwhm = ((img_table['fwhm_min_'+band][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
            center_coord = SkyCoord(img_table['RA_'+band][row], img_table['DEC_'+band][row], frame='icrs', unit=(u.deg, u.deg))
            center_coord_pix = center_coord.to_pixel(img_wcs)
            center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])

            pos_ang = (img_table['pa_'+band][row]-90)*u.deg
            
            ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm.value*2, pix_minor_fwhm.value*2, angle=pos_ang)
            size = pix_major_fwhm*2.1
            ap_mask = ellipse_reg.to_mask()
            cutout_mask = ap_mask.cutout(img_data)

            aperture_flux = np.nansum(cutout_mask[ap_mask.data==1])/ppbeam
            npix = len(cutout_mask[ap_mask.data==1])

            #now make annulus for measuring background and error
            annulus_width = 15 #pixels
            annulus_radius = img_table['fwhm_maj_'+band][row]*u.arcsecond#0.1*u.arcsecond
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
            pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
            bg_rms = median_abs_deviation(pixels_in_annulus)
            ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
            bg_median = np.median(pixels_in_annulus)

            pix_bg = bg_median*npix/ppbeam
            
            #flux corrections
            ap_flux_bgcorrect = aperture_flux - pix_bg
            ap_flux_correct = ap_flux_bgcorrect + ap_flux_bgcorrect*(1 - special.erf(2*np.sqrt(np.log(2)))) #flux correction for summing within 2*fwhm

            print(f'Background: {pix_bg}, Gauss amp: {img_table["gauss_amp_"+band][row]}')
            print(f'peak pixel: {np.nanmax(cutout_mask[ap_mask.data==1])}')

            
            ap_flux_err_arr.append(ap_bg_rms)
            ap_flux_arr.append(ap_flux_correct)
            snr_arr.append(img_table['gauss_amp_'+band][row]/bg_rms)

                
    cols = ['ap_flux_'+band, 'ap_flux_err_'+band, 'fwhm_maj_deconv_'+band, 'fwhm_maj_deconv_err_'+band, 'fwhm_min_deconv_'+band, 'fwhm_min_deconv_err_'+band, 'pa_deconv_'+band, 'pa_deconv_err_'+band, 'SNR_'+band]
    arrs = [ap_flux_arr, ap_flux_err_arr, fwhm_maj_deconv_arr, fwhm_maj_deconv_err_arr, fwhm_min_deconv_arr, fwhm_min_deconv_err_arr, pa_deconv_arr, pa_deconv_err_arr, snr_arr]
    for c in range(len(cols)):
        img_table.add_column(Column(np.array(arrs[c])), name=cols[c])
    img_table.add_column(Column(img_table['fwhm_maj_deconv_'+band]/img_table['fwhm_min_deconv_'+band]), name='ar_deconv_'+band)
    
    return img_table
