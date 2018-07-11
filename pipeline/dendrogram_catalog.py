#This file takes a dataset and dendrogram parameters and outputs a catalog of sources with gaussian fitted centroids
from region_file import compute_regions
from gaussfit_catalog import gaussfit_catalog
from regions import CircleSkyRegion

from astropy.io import fits, ascii
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, join
from astropy.nddata import Cutout2D

import numpy as np
import pylab as pl

#param_file = '/users/jotter/summer_research_2018/pipeline/data_params.csv'

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
		if params['pbcorr'][row] != 'True': #if the image is not pb corrected, measure flux from pb corrected image
			img = params['pbcorr'][row]
			corr_file = fits.open(img)
			data_corr = corr_file[0].data
			for i,leaf in enumerate(dendro.all_structures):
			    ind = leaf.indices()
			    max_y = np.max(ind[0])
			    min_y = np.min(ind[0])
			    max_x = np.max(ind[1])
			    min_x = np.min(ind[1])
			    cent = ((max_x + min_x)/2, (max_y + min_y)/2)
			    edges = (2*(max_x - min_x), 2*(max_y - min_y))
			    
			    cutout_corr = Cutout2D(data_corr,cent,edges)
			    
			    mask = leaf.get_mask()
			    mask_cut = Cutout2D(mask,cent,edges)
			    
			    leaf_data = cutout_corr.data[mask_cut.data==True]
			    leaf_flux = np.sum(leaf_data)
			    cat['flux'][leaf.idx] = leaf_flux
		cat.rename_column('flux', 'flux_'+name)
		cat.rename_column('_idx', '_idx_'+name)
		#next, measure centroids with gaussian fitting
		rad = Angle(1, 'arcsecond') #radius for region list
		#generate region list
		regs = []
		for ind in range(len(cat['x_cen'])):
		    regs.append(CircleSkyRegion(center=SkyCoord(cat['x_cen'][ind]*u.degree, cat['y_cen'][ind]*u.degree), radius=rad, meta={'text':str(cat['_idx_'+name][ind])}))

		cat_r = Angle(0.3, 'arcsecond') #radius in gaussian fitting
		pl.close('all')
		gauss_cat = gaussfit_catalog(img, regs, cat_r, savepath='/lustre/aoc/students/jotter/gauss_diags/'+name) #output is nested dictionary structure

		gauss_fit_tab = Table(names=('_idx_'+name, 'gauss_x_'+name, 'gauss_y_'+name, 'FWHM_major_'+name, 'FWHM_minor_'+name, 'position_angle_'+name, 'peak_flux_'+name, 'gauss_amplitude_'+name, 'fit_goodness_'+name), dtype=('i4','f8','f8','f8','f8','f8','f8','f8','U10')) #turn into astropy table
		for key in gauss_cat:
		    gauss_fit_tab.add_row((key,gauss_cat[key]['center_x'],gauss_cat[key]['center_y'], gauss_cat[key]['fwhm_major'], gauss_cat[key]['fwhm_minor'], gauss_cat[key]['pa'], gauss_cat[key]['peak'], gauss_cat[key]['amplitude'], 'none')) #fill table
	
		full_cat = join(cat, gauss_fit_tab, keys='_idx_'+name, join_type='outer') #joining the gaussian centroid data with the rest
	
		catalogs[name] = full_cat

	for dat in catalogs:
		Table(catalogs[dat]).write('/lustre/aoc/students/jotter/dendro_catalogs/'+dat+'_dendro_catalog_all.fits', format='fits', overwrite=True)

