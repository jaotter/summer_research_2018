#This file takes a dataset and dendrogram parameters and outputs a catalog of sources with gaussian fitted centroids
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

#param_file = '/users/jotter/summer_research_2018/pipeline/data_params.csv'

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

	catalogs = {}

	for row in range(len(params)):
		name=params['name'][row]
		min_val = params['min_value'][row]
		min_del = params['min_delta'][row]
		n_pix = params['n_pix'][row]
		img = params['file_name'][row]
		reg_fname = '/users/jotter/summer_research_2018/final_regs/'+name+'_reg_file.reg'
		dendro, cat = compute_regions(min_val, min_del, n_pix, img, reg_fname)
		cat['flux_err_'+name] = np.zeros(len(cat))
		if params['pbcorr'][row] != 'True': #if the image is not pb corrected, measure flux from pb corrected image
			img = params['pbcorr'][row]
			corr_file = fits.open(img)
			header = corr_file[0].header
			mywcs = WCS(header).celestial
			data = corr_file[0].data.squeeze()
			beam = radio_beam.Beam.from_fits_header(header)
			pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
			ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

			#for leaf in dendro.all_structures:
			leaf_inds = []
			for leaf in dendro.leaves:
			    leaf_inds.append(leaf.idx)
			    ind = leaf.indices()
			    max_y = np.max(ind[0])
			    min_y = np.min(ind[0])
			    max_x = np.max(ind[1])
			    min_x = np.min(ind[1])
			    cent = ((max_x + min_x)/2, (max_y + min_y)/2) #center of rectangular cutout
			    edges = (1.2*(max_x - min_x), 1.2*(max_y - min_y)) #1.2 chosen arbitrarily to give padding to cutout
			    
			    cutout_corr = Cutout2D(data,cent,edges)
			    
			    mask = leaf.get_mask()
			    mask_cut = Cutout2D(mask,cent,edges)
			    
			    leaf_data = cutout_corr.data[mask_cut.data==True]
			
			    npix = len(np.where(mask_cut.data == True)[0])
			    leaf_flux = np.sum(leaf_data)
			    flux_err = leaf_flux/np.sqrt(npix/ppbeam)
			    cat['flux'][leaf.idx] = leaf_flux
			    cat['flux_err_'+name][leaf.idx] = flux_err

		else:
			cont_file = fits.open(img)
			header = cont_file[0].header
			mywcs = WCS(header).celestial
			data = cont_file[0].data.squeeze()
			beam = radio_beam.Beam.from_fits_header(header)
			pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
			ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
			#for leaf in dendro.all_structures:
			leaf_inds = []
			for leaf in dendro.leaves:
				leaf_inds.append(leaf.idx)
				ind = leaf.indices()
				max_y = np.max(ind[0])
				min_y = np.min(ind[0])
				max_x = np.max(ind[1])
				min_x = np.min(ind[1])
				cent = ((max_x + min_x)/2, (max_y + min_y)/2) #center of rectangular cutout
				edges = (1.2*(max_x - min_x), 1.2*(max_y - min_y)) #1.2 chosen arbitrarily to give padding to cutout

				mask = leaf.get_mask()
				mask_cut = Cutout2D(mask,cent,edges)
				npix = len(np.where(mask_cut.data == True)[0])
				leaf_flux = cat['flux'][leaf.idx]
				flux_err = leaf_flux/np.sqrt(npix/ppbeam)
				cat['flux_err_'+name][leaf.idx] = flux_err
		cat.rename_column('flux', 'flux_'+name)
		cat.rename_column('_idx', '_idx_'+name)
		#next, measure centroids with gaussian fitting
		rad = Angle(1, 'arcsecond') #radius for region list
		#generate region list
		regs = []
		for ind in range(len(cat['x_cen'])):
		    if ind in leaf_inds:
		    	regs.append(regions.CircleSkyRegion(center=SkyCoord(cat['x_cen'][ind]*u.degree, cat['y_cen'][ind]*u.degree), radius=rad, meta={'text':str(cat['_idx_'+name][ind])}))

		cat_r = Angle(0.3, 'arcsecond') #radius in gaussian fitting
		#pl.close('all')
		gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath='/lustre/aoc/students/jotter/gauss_diags/leaves/'+name) #output is nested dictionary structure

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
			#cutout = Cutout2D(data, cutout_center, size.value, mywcs, mode='partial')
			ap_mask = ellipse_reg.to_mask()
			cutout_mask = ap_mask.cutout(data)

			aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
			gauss_fit_tab['ap_flux_'+name][row]=aperture_flux
			#npix = len(np.where(ap_mask.data == 1)[0])
			#ap_flux_err = aperture_flux/np.sqrt(npix/ppbeam)
			#gauss_fit_tab['ap_flux_err_'+name][row]=ap_flux_err
			

		full_cat = join(cat, gauss_fit_tab, keys='_idx_'+name, join_type='outer') #joining the gaussian centroid data with the rest
		full_cat.write('/lustre/aoc/students/jotter/dendro_catalogs/'+name+'_dendro_catalog_leaves.fits', format='fits', overwrite=True)


#data='/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'
#ref_name='B3'
#img = '/lustre/aoc/students/jotter/directory/Orion_SourceI_B7_continuum_r-2.mask5mJy.clean4mJy.image.tt0.pbcor.fits'
#name='B7_hr'

#img = '/lustre/aoc/students/jotter/directory/member.uid___A001_X88e_X1dd.Orion_BNKL_source_I_sci.spw25_27_29_31.cont.I.pbcor.fits'
#name='B7_lr'

def measure_fluxes(data, ref_name, img, name): #catalog to add measurements to, name of reference data set, img to get fluxes from, name of dataset
    cat = Table.read(data)

    fl = fits.open(img)
    header = fl[0].header
    img_data = fl[0].data.squeeze()

    mywcs = WCS(header).celestial

    rad = Angle(1, 'arcsecond') #radius for region list
    rad_deg = rad.to(u.degree)
    #generate region list
    regs = []

    reg_file = '/users/jotter/summer_research_2018/final_regs/'+name+'_reg_file.reg'

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
    #pl.close('all')
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
        ellipse_reg = regions.EllipsePixelRegion(center_coord_pix, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
        size = pix_major_fwhm*2.1
        ap_mask = ellipse_reg.to_mask()
        cutout_mask = ap_mask.cutout(img_data)

        aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
        gauss_fit_tab['ap_flux_'+name][row]=aperture_flux
    
    gauss_fit_tab.write('/lustre/aoc/students/jotter/dendro_catalogs/'+name+'_dendro_catalog_ref_'+ref_name+'.fits', format='fits', overwrite=True)


#data = '/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits'
#directory = '/lustre/aoc/students/jotter/directory/OrionB6/'
#ref_name = 'B6'
#name_start = 27
#name_end = -21
def multi_cat_measure_fluxes(data, ref_name, directory, name_start = 27, name_end = -21): #catalog to get positions from, name of reference data, directory of images to get fluxes from 
    cat = Table.read(data)
    imgs = glob.glob(directory+'*')
    img_names = [img[len(directory)+name_start:name_end] for img in imgs]

    for j,img in enumerate(imgs):
        name = img_names[j]
        if len(name) > 54:
            name = name[0:53]
        fl = fits.open(img)
        header = fl[0].header
        img_data = fl[0].data.squeeze()

        mywcs = WCS(header).celestial

        rad = Angle(1, 'arcsecond') #radius for region list
        rad_deg = rad.to(u.degree)
        #generate region list
        regs = []
        
        save_dir = '/lustre/aoc/students/jotter/gauss_diags/diff_imgs/'+ref_name+'/'+name+'/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        for ind in range(len(cat)):
            reg = regions.CircleSkyRegion(center=SkyCoord(cat['gauss_x_'+ref_name][ind]*u.degree, cat['gauss_y_'+ref_name][ind]*u.degree), radius=rad, meta={'text':str(cat['D_ID'][ind])})
            reg_pix = reg.to_pixel(mywcs)
            if reg_pix.center.x > 0 and reg_pix.center.x < len(img_data[0]):
                if reg_pix.center.y > 0 and reg_pix.center.y < len(img_data):
                    regs.append(reg)

        cat_r = Angle(0.5, 'arcsecond') #radius in gaussian fitting
        gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath=save_dir) #output is nested dictionary structure

        gauss_fit_tab = Table(names=('D_ID','FWHM_major_'+name, 'major_err_'+name, 'FWHM_minor_'+name, 'minor_err_'+name, 'pa_'+name, 'pa_err_'+name, 'g_amplitude_'+name, 'amp_err_'+name, 'ap_flux_'+name, 'ap_flux_err_'+name, 'bg_median_'+name), dtype=('i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')) 
        for key in gauss_cat:
            gauss_fit_tab.add_row((key, gauss_cat[key]['fwhm_major'], gauss_cat[key]['e_fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['e_fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['e_pa'], gauss_cat[key]['amplitude'], gauss_cat[key]['e_amplitude'], np.nan, np.nan, np.nan)) #fill table

        beam = radio_beam.Beam.from_fits_header(header)
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
        ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

        for row in range(len(gauss_fit_tab)):
            pix_major_fwhm = ((gauss_fit_tab['FWHM_major_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()
            pix_minor_fwhm = ((gauss_fit_tab['FWHM_minor_'+name][row]*u.arcsec).to(u.degree)/pixel_scale).decompose()

            cat_ind = np.where(cat['D_ID'] == gauss_fit_tab['D_ID'][row])
            center_coord = SkyCoord(cat['gauss_x_'+ref_name][cat_ind], cat['gauss_y_'+ref_name][cat_ind], frame='icrs', unit=(u.deg, u.deg))
            center_coord_pix = center_coord.to_pixel(mywcs)
            center_coord_pix_reg = regions.PixCoord(center_coord_pix[0], center_coord_pix[1])
            position_angle = gauss_fit_tab['pa_'+name][row]*u.deg
            ellipse_reg = regions.EllipsePixelRegion(center_coord_pix_reg, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
            size = pix_major_fwhm*2.1
            ap_mask = ellipse_reg.to_mask()
            cutout_mask = ap_mask.cutout(img_data)

            aperture_flux = np.sum(cutout_mask[ap_mask.data==1])/ppbeam
            gauss_fit_tab['ap_flux_'+name][row]=aperture_flux

            #now creating annulus around source to measure background and ap flux error
            annulus_width = 15
            center_pad = 10 #pad between ellipse and inner radius

            # Cutout section of the image we care about, to speed up computation time
            size = ((center_pad+annulus_width)+pix_major_fwhm)*2.2*pixel_scale #2.2 is arbitrary to get entire annulus and a little extra
            cutout = Cutout2D(img_data, center_coord_pix, size, mywcs, mode='partial') #cutout of outer circle
            cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])


            # Define the aperture regions needed for SNR
            innerann_reg = regions.CirclePixelRegion(cutout_center, center_pad+pix_major_fwhm)
            outerann_reg = regions.CirclePixelRegion(cutout_center, center_pad+pix_major_fwhm+annulus_width)

            # Make masks from aperture regions
            annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)
            
            # Calculate the SNR and aperture flux sums
            pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
            bg_rms = rms(pixels_in_annulus)
            bg_median = np.median(pixels_in_annulus)

            gauss_fit_tab['ap_flux_err_'+name][row] = bg_rms
            gauss_fit_tab['bg_median_'+name][row] = bg_median
        cat = join(cat, gauss_fit_tab, keys='D_ID', join_type='left')
        
    cat.write('/lustre/aoc/students/jotter/dendro_catalogs/'+ref_name+'_diff_imgs_catalog.txt', format='ascii', overwrite=True)

