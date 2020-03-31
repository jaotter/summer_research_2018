from region_file import compute_regions
from gaussfit_catalog import gaussfit_catalog
import regions
import radio_beam

from astropy.io import fits, ascii
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, join
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

import numpy as np
import pylab as pl
import glob
import os
import fnmatch
from functools import reduce

#This file takes a dataset and dendrogram parameters and outputs a catalog of sources with gaussian fitted centroids
#Dependencies: 'region_file.py' in the same directory. There are also paths to where region files and gaussian fit diagnostic plots should be saved, which may need to be changed.

#create_catalog is the main function, and it requires a file of parameters as input. For an example, see 'param_file.csv'. 

#The other two main functions, measure_ap_fluxes and measure_fluxes_gaussfit are meant for entire directories of image files
#measure_ap_fluxes adds the aperture fluxes of all files in a directory to a single catalog
#measure_fluxes_gaussfit does the same as above but also fits gaussians, which takes much longer

def rms(array):
	sq_arr = np.square(array)
	avg = np.nanmean(sq_arr)
	return np.sqrt(avg)

def mask(reg, cutout):#masks everything except the region
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int') 


def create_catalog(param_file):
	params = ascii.read(param_file)

	for row in range(len(params)):
		name=params['name'][row] #contains parameters for dendrogram catalog creation
		min_val = params['min_value'][row]
		min_del = params['min_delta'][row]
		n_pix = params['n_pix'][row]
		img = params['file_name'][row]
		reg_fname = '/home/jotter/nrao/images/dendro_regions/'+name+'_reg_file.reg'
		dendro, cat = compute_regions(min_val, min_del, n_pix, img, reg_fname)
		cat['flux_err_'+name] = np.zeros(len(cat))
		if params['pbcorr'][row] != 'True': #if the image is not pb corrected, measure flux from pb corrected image
			img = params['pbcorr'][row]
		cont_file = fits.open(img)
		header = cont_file[0].header
		mywcs = WCS(header).celestial
		data = cont_file[0].data.squeeze()
		beam = radio_beam.Beam.from_fits_header(header)
		pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
		ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

		cat.rename_column('_idx', '_idx_'+name)

		#next, measure centroids with gaussian fitting
		rad = Angle(1, 'arcsecond') #radius for region list
		
		regs = []
		for ind,leaf in enumerate(dendro.all_structures): #.leaves #generate region list
			regs.append(regions.CircleSkyRegion(center=SkyCoord(cat['x_cen'][ind]*u.degree, cat['y_cen'][ind]*u.degree), radius=rad, meta={'text':str(cat['_idx_'+name][ind])}))
		print(len(regs))
		#if name == 'B6':
                        #regs.append(regions.CircleSkyRegion(center=SkyCoord(83.81138026377482*u.degree, -5.374951161716349*u.degree), radius=rad, meta={'text':'54'}))#manually appending a source not picked up by dendrogram

		cat_r = Angle(0.3, 'arcsecond') #radius in gaussian fitting
		gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath='/home/jotter/nrao/gauss_diags/leaves/'+name) #output is nested dictionary structure

		gauss_fit_tab = Table(names=('_idx_'+name, 'gauss_x_'+name, 'x_err_'+name, 'gauss_y_'+name, 'y_err_'+name, 'FWHM_major_'+name, 'major_err_'+name, 'FWHM_minor_'+name, 'minor_err_'+name, 'position_angle_'+name, 'position_angle_err_'+name, 'peak_flux_'+name, 'gauss_amplitude_'+name, 'amplitude_err_'+name, 'ap_flux_'+name, 'ap_flux_err_'+name,'fit_goodness_'+name), dtype=('i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','U10')) #turn into astropy table
		for key in gauss_cat:
		    gauss_fit_tab.add_row((key,gauss_cat[key]['center_x'],gauss_cat[key]['e_center_x'],gauss_cat[key]['center_y'],gauss_cat[key]['e_center_y'], gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['peak'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], np.nan, np.nan,'none')) #fill table
	
		for row in range(len(gauss_fit_tab)):
			pix_major_fwhm = ((gauss_fit_tab['FWHM_major_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
			pix_minor_fwhm = ((gauss_fit_tab['FWHM_minor_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
			cutout_center = SkyCoord(gauss_fit_tab['gauss_x_'+name][row], gauss_fit_tab['gauss_y_'+name][row], frame='icrs', unit=(u.deg, u.deg))
			cutout_center_pix = cutout_center.to_pixel(mywcs)
			cutout_center_pix = regions.PixCoord(cutout_center_pix[0], cutout_center_pix[1])
			position_angle = gauss_fit_tab['position_angle_'+name][row]*u.deg
			ellipse_reg = regions.EllipsePixelRegion(cutout_center_pix, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
			size = pix_major_fwhm*2.1
			ap_mask = ellipse_reg.to_mask()
			cutout_mask = ap_mask.cutout(data)

			aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
			gauss_fit_tab['ap_flux_'+name][row]=aperture_flux

			#NOTE: aperture flux error and background fluxes will get filled in after running `snr_rejection.py`

		full_cat = join(cat, gauss_fit_tab, keys='_idx_'+name, join_type='outer') #joining the gaussian centroid data with the rest
		full_cat.write('/home/jotter/nrao/tables/dendro_catalogs/'+name+'_dendro_catalog_leaves.fits', format='fits', overwrite=True)


def measure_fluxes(data, ref_name, img, name):

    #This function is to add another image to an already made catalog using positions from a different image. 'data' is the catalog to add measurements to, 'ref_name' is the name of the reference data set (where the positions are from), 'img' is the new image, and 'name' is the name of 'img'

    cat = Table.read(data)

    fl = fits.open(img)
    header = fl[0].header
    img_data = fl[0].data.squeeze()

    mywcs = WCS(header).celestial

    rad = Angle(1, 'arcsecond') #radius for region list
    rad_deg = rad.to(u.degree)
    #generate region list
    regs = []

    reg_file = '/users/jotter/summer_research_2018/final_regs/'+name+'_reg_file_circle.reg'

    with open(reg_file, 'w') as fh:
        fh.write("fk5\n")
        for ind in range(len(cat)):
            reg = regions.CircleSkyRegion(center=SkyCoord(cat['gauss_x_'+ref_name][ind]*u.degree, cat['gauss_y_'+ref_name][ind]*u.degree), radius=rad, meta={'text':str(cat['D_ID'][ind])})
            reg_pix = reg.to_pixel(mywcs)
            if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                    regs.append(reg)
                    fh.write('circle({x_cen}, {y_cen}, {radius}) #text={{{ID}}}\n'.format(x_cen=reg.center.ra.value, y_cen=reg.center.dec.value, radius=0.1*rad_deg.value, ID=str(cat['D_ID'][ind])))

    cat_r = Angle(0.3, 'arcsecond') #radius in gaussian fitting
    gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath='/lustre/aoc/students/jotter/gauss_diags/leaves/'+name) #output is nested dictionary structure

    gauss_fit_tab = Table(names=('D_ID', '_idx_'+name, 'gauss_x_'+name, 'x_err_'+name, 'gauss_y_'+name, 'y_err_'+name, 'FWHM_major_'+name, 'major_err_'+name, 'FWHM_minor_'+name, 'minor_err_'+name, 'position_angle_'+name, 'position_angle_err_'+name, 'peak_flux_'+name, 'gauss_amplitude_'+name, 'amplitude_err_'+name, 'ap_flux_'+name, 'ap_flux_err_'+name,'fit_goodness_'+name), dtype=('i4','i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','U10')) #turn into astropy table
    for key in gauss_cat:
        gauss_fit_tab.add_row((key, key, gauss_cat[key]['center_x'],gauss_cat[key]['e_center_x'],gauss_cat[key]['center_y'],gauss_cat[key]['e_center_y'], gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['peak'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], np.nan, np.nan,'none')) #fill table

    beam = radio_beam.Beam.from_fits_header(header)
    pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
    ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

    for row in range(len(gauss_fit_tab)):
        pix_major_fwhm = ((gauss_fit_tab['FWHM_major_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
        pix_minor_fwhm = ((gauss_fit_tab['FWHM_minor_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
        cutout_center = SkyCoord(gauss_fit_tab['gauss_x_'+name][row], gauss_fit_tab['gauss_y_'+name][row], frame='icrs', unit=(u.deg, u.deg))
        cutout_center_pix = cutout_center.to_pixel(mywcs)
        cutout_center_pix = regions.PixCoord(cutout_center_pix[0], cutout_center_pix[1])
        position_angle = gauss_fit_tab['position_angle_'+name][row]*u.deg
        ellipse_reg = regions.EllipsePixelRegion(cutout_center_pix, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
        size = pix_major_fwhm*2.1
        ap_mask = ellipse_reg.to_mask()
        cutout_mask = ap_mask.cutout(img_data)

        aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
        gauss_fit_tab['ap_flux_'+name][row]=aperture_flux
    
    gauss_fit_tab.write('/lustre/aoc/students/jotter/dendro_catalogs/'+name+'_dendro_catalog_ref_'+ref_name+'.fits', format='fits', overwrite=True)


def measure_ap_fluxes(data, directory, data_name, name_start = 27, name_end = -21, ref_name = 'B6'): 

    #This function takes a directory full of images and adds them all to a catalog with flux measurements. 'data' is the catalog to get positions and sizes from, 'ref_name' is the name of the reference data (for positions), 'directory' is the directory with all the new images.
    #data = '/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits'
    #directory = '/lustre/aoc/students/jotter/directory/OrionB3/'

    cat = Table.read(data)
    RA_names = ['gauss_x_B3', 'gauss_x_B6', 'gauss_x_B7_hr']
    inds = []
    for fn in RA_names:
        ind = np.where(np.isnan(cat[fn]) == False)
        inds.append(ind)
    detected = reduce(np.intersect1d, inds) #sources detected in B3, B6, B7 - only make measurements on these sources
    cat = cat[detected]

    imgs = glob.glob(directory+'*')
    img_names = [img[len(directory)+name_start:name_end] for img in imgs]

    col_names = fnmatch.filter(cat.colnames, 'ap_flux_r*')
    col_names = [cn[8:] for cn in col_names]

    for j,img in enumerate(imgs):
        name = img_names[j]
        if len(name) > 54:
            name = name[0:53]
        if name not in col_names:
            fl = fits.open(img)
            header = fl[0].header
            try:
                img_data = fl[0].data.squeeze()
            except TypeError:
                print("ERROR! img: "+img)
                continue
            mywcs = WCS(header).celestial

            beam = radio_beam.Beam.from_fits_header(header)
            pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
            ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

            ap_flux_arr = []
            circ_flux_arr = []
            ap_flux_err_arr = []
            circ_flux_err_arr = []
            bg_median_arr = []
            bg_ap_arr = []
            bg_circ_arr = []
            for row in range(len(cat)):
                pix_major_fwhm = ((cat['FWHM_major_'+ref_name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                pix_minor_fwhm = ((cat['FWHM_minor_'+ref_name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()

                center_coord = SkyCoord(cat['gauss_x_'+ref_name][row], cat['gauss_y_'+ref_name][row], frame='icrs', unit=(u.deg, u.deg))
                center_coord_pix = center_coord.to_pixel(mywcs)
                center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])
                position_angle = cat['position_angle_'+ref_name][row]*u.deg
                ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
                size = pix_major_fwhm*2.1
                ap_mask = ellipse_reg.to_mask()
                cutout_mask = ap_mask.cutout(img_data)

                aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
                npix = len(cutout_mask[ap_mask.data==1])
                ap_flux_arr.append(aperture_flux)

                #now creating annulus around source to measure background and ap flux error
                annulus_width = 15
                annulus_radius = 0.1*u.arcsecond
                annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

                # Cutout section of the image we care about, to speed up computation time
                size = 2.5*annulus_radius
                cutout = Cutout2D(img_data, center_coord_pix, size, mywcs, mode='partial') #cutout of outer circle
                cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

                # Define the aperture regions needed for SNR
                innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix)
                outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix+annulus_width)

                # Make masks from aperture regions
                annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)
                
                # Calculate the SNR and aperture flux sums
                pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
                bg_rms = rms(pixels_in_annulus)
                ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
                bg_median = np.median(pixels_in_annulus)

                pix_bg = bg_median*npix/ppbeam

                ap_flux_err_arr.append(ap_bg_rms)
                bg_median_arr.append(bg_median)
                bg_ap_arr.append(pix_bg)

                #now measure circle flux:
                radius = 0.1*u.arcsecond
                radius_pix = annulus_radius_pix
                circle_reg = regions.CirclePixelRegion(center_coord_pix_reg, radius_pix)
                circ_ap_mask = circle_reg.to_mask()
                circ_cutout_mask = circ_ap_mask.cutout(img_data)
                cutout_mask = ap_mask.cutout(img_data)
                circ_aperture_flux = np.sum(circ_cutout_mask[circ_ap_mask.data==1])/ppbeam
                circ_npix = len(circ_cutout_mask[circ_ap_mask.data==1])

                circ_bg_rms = bg_rms/np.sqrt(circ_npix/ppbeam)

                circ_flux_arr.append(circ_aperture_flux)
                circ_flux_err_arr.append(circ_bg_rms)
                bg_circ_arr.append(bg_median*circ_npix/ppbeam)
                

            cols = ['ap_flux_', 'ap_flux_err_', 'bg_median_', 'bg_ap_', 'circ_flux_', 'circ_flux_err_', 'bg_circ_']
            arrs = [ap_flux_arr, ap_flux_err_arr, bg_median_arr, bg_ap_arr, circ_flux_arr, circ_flux_err_arr, bg_circ_arr]
            for j in range(len(cols)):
                cat[cols[j]+name] = arrs[j]

        
    cat.write('/lustre/aoc/students/jotter/dendro_catalogs/'+data_name+'_diff_imgs_catalog_'+ref_name+'_ref.txt', format='ascii', overwrite=True)



def measure_fluxes_gaussfit(data, directory, data_name, name_start = 27, name_end = -21, ref_name = 'B6'): 

    #This function takes a directory full of images and adds them all to a catalog with flux measurements and fits gaussians. 'data' is the catalog to get positions from, 'ref_name' is the name of the reference data (for positions), 'directory' is the directory with all the new images.

    cat = Table.read(data)
    RA_names = ['gauss_x_B3', 'gauss_x_B6', 'gauss_x_B7_hr']
    inds = []
    for fn in RA_names:
        ind = np.where(np.isnan(cat[fn]) == False)
        inds.append(ind)
    detected = reduce(np.intersect1d, inds) #sources detected in B3, B6, B7 - only make measurements on these sources
    cat = cat[detected]

    imgs = glob.glob(directory+'*r-2.clean0.1mJy.500klplus.deepmask*')
    img_names = [img[len(directory)+name_start:name_end] for img in imgs]

    col_names = fnmatch.filter(cat.colnames, 'ap_flux_r*')
    col_names = [cn[8:] for cn in col_names]

    for j,img in enumerate(imgs):
        name = img_names[j]
        if len(name) > 54:
            name = name[0:53]
        if name not in col_names:
            fl = fits.open(img)
            header = fl[0].header
            try:
                img_data = fl[0].data.squeeze()
            except TypeError:
                print("error! img: "+img)
                continue
            mywcs = WCS(header).celestial

            beam = radio_beam.Beam.from_fits_header(header)
            pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
            ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

            save_dir = '/lustre/aoc/students/jotter/gauss_diags/diff_imgs/'+ref_name+'/'+name+'/'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
    
            reg_file = '/users/jotter/summer_research_2018/final_regs/'+name+'_reg_file_B6_apflux.reg'
            with open(reg_file, 'w') as fh:
                fh.write("fk5\n")
                for ind in range(len(cat)):
                    fh.write('ellipse({x_cen}, {y_cen}, {maj}, {minr}, {ang}) #text={{{ID}}}\n'.format(x_cen=cat['gauss_x_'+ref_name][ind], y_cen=cat['gauss_y_'+ref_name][ind], maj = (cat['FWHM_major_'+ref_name][ind]*u.arcsec.to(u.degree)), minr = (cat['FWHM_minor_'+ref_name][ind]*u.arcsec.to(u.degree)), ang = cat['position_angle_'+ref_name][ind], ID=str(cat['D_ID'][ind])))

            rad = Angle(1, 'arcsecond') #radius for region list
            regs = []

            for ind in range(len(cat)):
                reg = regions.CircleSkyRegion(center=SkyCoord(cat['gauss_x_'+ref_name][ind]*u.degree, cat['gauss_y_'+ref_name][ind]*u.degree), radius=rad, meta={'text':str(cat['D_ID'][ind])})
                reg_pix = reg.to_pixel(mywcs)
                if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                    if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                        regs.append(reg)

            cat_r = Angle(0.5, 'arcsecond') #radius in gaussian fitting
            gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=save_dir) #output is nested dictionary structure

            gauss_fit_tab = Table(names=('D_ID','FWHM_major_'+name, 'major_err_'+name, 'FWHM_minor_'+name, 'minor_err_'+name, 'pa_'+name, 'pa_err_'+name, 'g_amplitude_'+name, 'amp_err_'+name), dtype=('i4','f8','f8','f8','f8','f8','f8','f8','f8')) 
            for key in gauss_cat:
                gauss_fit_tab.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'])) #fill table 

            cat = join(cat, gauss_fit_tab, keys='D_ID', join_type='left')
            
            ap_flux_arr = []
            circ_flux_arr = []
            ap_flux_err_arr = []
            circ_flux_err_arr = []
            bg_median_arr = []
            bg_ap_arr = []
            bg_circ_arr = []

            for row in range(len(cat)):
                pix_major_fwhm = ((cat['FWHM_major_'+ref_name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                pix_minor_fwhm = ((cat['FWHM_minor_'+ref_name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()

                center_coord = SkyCoord(cat['gauss_x_'+ref_name][row], cat['gauss_y_'+ref_name][row], frame='icrs', unit=(u.deg, u.deg))
                center_coord_pix = center_coord.to_pixel(mywcs)
                center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])
                position_angle = cat['position_angle_'+ref_name][row]*u.deg
                print(center_coord_pix_reg, pix_major_fwhm, pix_minor_fwhm)
                ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
                size = pix_major_fwhm*2.1
                ap_mask = ellipse_reg.to_mask()
                cutout_mask = ap_mask.cutout(img_data)

                aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
                npix = len(cutout_mask[ap_mask.data==1])
                ap_flux_arr.append(aperture_flux)

                #now creating annulus around source to measure background and ap flux error
                annulus_width = 15
                annulus_radius = 0.1*u.arcsecond
                annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

                # Cutout section of the image we care about, to speed up computation time
                size = 2.5*annulus_radius
                cutout = Cutout2D(img_data, center_coord_pix, size, mywcs, mode='partial') #cutout of outer circle
                cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

                # Define the aperture regions needed for SNR
                innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix)
                outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix+annulus_width)

                # Make masks from aperture regions
                annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)
                
                # Calculate the SNR and aperture flux sums
                pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
                bg_rms = rms(pixels_in_annulus)
                ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
                bg_median = np.median(pixels_in_annulus)

                pix_bg = bg_median*npix/ppbeam

                ap_flux_err_arr.append(ap_bg_rms)
                bg_median_arr.append(bg_median)
                bg_ap_arr.append(pix_bg)

                #now measure circle flux:
                radius = 0.1*u.arcsecond
                radius_pix = annulus_radius_pix
                circle_reg = regions.CirclePixelRegion(center_coord_pix_reg, radius_pix)
                circ_ap_mask = circle_reg.to_mask()
                circ_cutout_mask = circ_ap_mask.cutout(img_data)
                cutout_mask = ap_mask.cutout(img_data)
                circ_aperture_flux = np.sum(circ_cutout_mask[circ_ap_mask.data==1])/ppbeam
                circ_npix = len(circ_cutout_mask[circ_ap_mask.data==1])

                circ_bg_rms = bg_rms/np.sqrt(circ_npix/ppbeam)

                circ_flux_arr.append(circ_aperture_flux)
                circ_flux_err_arr.append(circ_bg_rms)
                bg_circ_arr.append(bg_median*circ_npix/ppbeam)
                

            cols = ['ap_flux_', 'ap_flux_err_', 'bg_median_', 'bg_ap_', 'circ_flux_', 'circ_flux_err_', 'bg_circ_']
            arrs = [ap_flux_arr, ap_flux_err_arr, bg_median_arr, bg_ap_arr, circ_flux_arr, circ_flux_err_arr, bg_circ_arr]
            for j in range(len(cols)):
                cat[cols[j]+name] = arrs[j]
        
    cat.write('/lustre/aoc/students/jotter/dendro_catalogs/'+data_name+'_500klplus_allsrcs_catalog_'+ref_name+'_ref.txt', format='ascii', overwrite=True)


def gauss_sizes_imgs(srcID, band, images, directory, name_start = 27, name_end = -21): 

    #This function takes a list of image names and adds them all to a catalog with gaussian sizes. 'data' is the catalog to get positions from, 'srcID' is the D_ID of the desired source, 'band' is the band of data, 'images' is the list of image names, 'directory' is the directory where these images are.

    data = '/lustre/aoc/students/jotter/dendro_catalogs/master_500klplus_B3_ref.fits'
    cat = Table.read(data)
    ind = np.where(cat['D_ID'] == srcID)[0][0]
    
    img_names = [img[name_start:name_end] for img in images]
    print(img_names)
    
    table_names = [] #because maybe not all images in img_namse get used
    img_fwhm_maj = [] #at the end combine these lists into a table
    img_fwhm_maj_err = []
    img_fwhm_min = []
    img_fwhm_min_err = []
    img_pa = []
    img_pa_err = []
    
    for j,img in enumerate(images): #loop through images and measure fwhms
        name = img_names[j]
        if len(name) > 54:
            name = name[0:53]

        fl = fits.open(directory+img)
        header = fl[0].header
        try:
            img_data = fl[0].data.squeeze()
        except TypeError:
            print("error! img: "+img)
            continue
        mywcs = WCS(header).celestial

        beam = radio_beam.Beam.from_fits_header(header)
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
        ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

        save_dir = '/lustre/aoc/students/jotter/gauss_diags/diff_imgs/src_'+str(srcID)+'/'+band+'/'+name
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        rad = Angle(1, 'arcsecond') #radius for region list
        regs = []

        reg = regions.CircleSkyRegion(center=SkyCoord(cat['RA_B3'][ind]*u.degree, cat['DEC_B3'][ind]*u.degree), radius=rad, meta={'text':str(cat['D_ID'][ind])})
        reg_pix = reg.to_pixel(mywcs)
        if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
            if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                regs.append(reg)

        cat_r = Angle(0.5, 'arcsecond') #radius in gaussian fitting
        gauss_cat = gaussfit_catalog(directory+img, regs, cat_r, savepath=save_dir) #output is nested dictionary structure

        table_names.append(name)
        img_fwhm_maj.append(gauss_cat[str(srcID)]['fwhm_major'].value)
        img_fwhm_maj_err.append(gauss_cat[str(srcID)]['e_fwhm_major'].value)
        img_fwhm_min.append(gauss_cat[str(srcID)]['fwhm_minor'].value)
        img_fwhm_min_err.append(gauss_cat[str(srcID)]['e_fwhm_minor'].value)
        img_pa.append(gauss_cat[str(srcID)]['pa'].value)
        img_pa_err.append(gauss_cat[str(srcID)]['e_pa'].value)

    table = Table((table_names, img_fwhm_maj, img_fwhm_maj_err, img_fwhm_min, img_fwhm_min_err, img_pa, img_pa_err), names=('img_name', 'fwhm_maj', 'fwhm_maj_err', 'fwhm_min', 'fwhm_min_err', 'pa', 'pa_err'))
    table['fwhm_maj'].unit = 'arcsec'
    table['fwhm_maj_err'].unit = 'arcsec'
    table['fwhm_min'].unit = 'arcsec'
    table['fwhm_min_err'].unit = 'arcsec'

    table.write('/lustre/aoc/students/jotter/dendro_catalogs/src'+str(srcID)+'_img_sizes'+band+'.fits', overwrite=True)

def measure_fluxes_gaussfit_allsrcs(data, img, data_name, name_start = 27, name_end = -21, ref_name = 'B3'): 

    cat = Table.read(data)
    RA_name = 'gauss_x_'+ref_name

    ind = np.where(np.isnan(cat[RA_name]) == False)
    cat = cat[ind]

    col_names = fnmatch.filter(cat.colnames, 'ap_flux_r*')
    col_names = [cn[8:] for cn in col_names]

    fl = fits.open(img)
    header = fl[0].header
    img_data = fl[0].data.squeeze()
    mywcs = WCS(header).celestial

    beam = radio_beam.Beam.from_fits_header(header)
    pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
    ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

    save_dir = '/lustre/aoc/students/jotter/gauss_diags/diff_imgs/'+ref_name+'/'+data_name+'/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    reg_file = '/users/jotter/summer_research_2018/final_regs/misc_regs/'+data_name+'_reg_file_'+ref_name+'_apflux.reg'
    with open(reg_file, 'w') as fh:
        fh.write("fk5\n")
        for ind in range(len(cat)):
            fh.write('ellipse({x_cen}, {y_cen}, {maj}, {minr}, {ang}) #text={{{ID}}}\n'.format(x_cen=cat['gauss_x_'+ref_name][ind], y_cen=cat['gauss_y_'+ref_name][ind], maj = (cat['FWHM_major_'+ref_name][ind]*u.arcsec.to(u.degree)), minr = (cat['FWHM_minor_'+ref_name][ind]*u.arcsec.to(u.degree)), ang = cat['position_angle_'+ref_name][ind], ID=str(cat['D_ID'][ind])))

    rad = Angle(1, 'arcsecond') #radius for region list
    regs = []

    for ind in range(len(cat)):
        reg = regions.CircleSkyRegion(center=SkyCoord(cat['gauss_x_'+ref_name][ind]*u.degree, cat['gauss_y_'+ref_name][ind]*u.degree), radius=rad, meta={'text':str(cat['D_ID'][ind])})
        reg_pix = reg.to_pixel(mywcs)
        if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
            if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                if np.isnan(img_data[int(reg_pix.center.x), int(reg_pix.center.y)]) == False:
                    regs.append(reg)

    cat_r = Angle(0.5, 'arcsecond') #radius in gaussian fitting
    gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=save_dir) #output is nested dictionary structure

    gauss_fit_tab = Table(names=('D_ID','FWHM_major_'+data_name, 'major_err_'+data_name, 'FWHM_minor_'+data_name, 'minor_err_'+data_name, 'pa_'+data_name, 'pa_err_'+data_name, 'g_amplitude_'+data_name, 'amp_err_'+data_name), dtype=('i4','f8','f8','f8','f8','f8','f8','f8','f8')) 
    for key in gauss_cat:
        gauss_fit_tab.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'])) #fill table 
    
    ap_flux_arr = []
    circ_flux_arr = []
    ap_flux_err_arr = []
    circ_flux_err_arr = []
    bg_median_arr = []
    bg_ap_arr = []
    bg_circ_arr = []
    RA = []
    DEC = []
    RA_err = []
    DEC_err = []

    for row in range(len(cat)):
        gauss_ind = np.where(gauss_fit_tab['D_ID'] == cat['D_ID'][row])[0]
        if len(gauss_ind) > 0:
            pix_major_fwhm = ((cat['FWHM_major_'+ref_name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
            pix_minor_fwhm = ((cat['FWHM_minor_'+ref_name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()

            center_coord = SkyCoord(cat['gauss_x_'+ref_name][row], cat['gauss_y_'+ref_name][row], frame='icrs', unit=(u.deg, u.deg))
            center_coord_pix = center_coord.to_pixel(mywcs)
            center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])
            position_angle = cat['position_angle_'+ref_name][row]*u.deg

            ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
            size = pix_major_fwhm*2.1
            ap_mask = ellipse_reg.to_mask()
            cutout_mask = ap_mask.cutout(img_data)

            aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
            npix = len(cutout_mask[ap_mask.data==1])
            ap_flux_arr.append(aperture_flux)

            #now creating annulus around source to measure background and ap flux error
            annulus_width = 15
            annulus_radius = 0.1*u.arcsecond
            annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

            # Cutout section of the image we care about, to speed up computation time
            size = 2.5*annulus_radius
            cutout = Cutout2D(img_data, center_coord_pix, size, mywcs, mode='partial') #cutout of outer circle
            cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

            # Define the aperture regions needed for SNR
            innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix)
            outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix+annulus_width)

            # Make masks from aperture regions
            annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)
            
            # Calculate the SNR and aperture flux sums
            pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
            bg_rms = rms(pixels_in_annulus)
            ap_bg_rms = bg_rms/np.sqrt(npix/ppbeam) #rms/sqrt(npix/ppbeam) - rms error per beam
            bg_median = np.median(pixels_in_annulus)

            pix_bg = bg_median*npix/ppbeam

            ap_flux_err_arr.append(ap_bg_rms)
            bg_median_arr.append(bg_median)
            bg_ap_arr.append(pix_bg)

            #now measure circle flux:
            radius = 0.1*u.arcsecond
            radius_pix = annulus_radius_pix
            circle_reg = regions.CirclePixelRegion(center_coord_pix_reg, radius_pix)
            circ_ap_mask = circle_reg.to_mask()
            circ_cutout_mask = circ_ap_mask.cutout(img_data)
            cutout_mask = ap_mask.cutout(img_data)
            circ_aperture_flux = np.sum(circ_cutout_mask[circ_ap_mask.data==1])/ppbeam
            circ_npix = len(circ_cutout_mask[circ_ap_mask.data==1])

            circ_bg_rms = bg_rms/np.sqrt(circ_npix/ppbeam)

            circ_flux_arr.append(circ_aperture_flux)
            circ_flux_err_arr.append(circ_bg_rms)
            bg_circ_arr.append(bg_median*circ_npix/ppbeam)

            RA.append(cat['gauss_x_'+ref_name][row])
            DEC.append(cat['gauss_y_'+ref_name][row])
            RA_err.append(cat['x_err_'+ref_name][row])
            DEC_err.append(cat['y_err_'+ref_name][row])
        

    cols = ['ap_flux_', 'ap_flux_err_', 'bg_median_', 'bg_ap_', 'circ_flux_', 'circ_flux_err_', 'bg_circ_']
    arrs = [ap_flux_arr, ap_flux_err_arr, bg_median_arr, bg_ap_arr, circ_flux_arr, circ_flux_err_arr, bg_circ_arr]
    cols2 = ['RA_', 'DEC_', 'RA_err_', 'DEC_err_'] #seperate bc different naming
    arrs2 = [RA, DEC, RA_err, DEC_err]
    for j in range(len(cols)):
        gauss_fit_tab[cols[j]+data_name] = arrs[j]
    for c in range(len(cols2)):
        gauss_fit_tab[cols2[c]+ref_name] = arrs2[c]
        
    gauss_fit_tab.write('/lustre/aoc/students/jotter/dendro_catalogs/'+data_name+'_500klplus_allsrcs_catalog_'+ref_name+'_ref.txt', format='ascii', overwrite=True)

