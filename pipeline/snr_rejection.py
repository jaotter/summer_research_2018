from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
import numpy as np
import regions
from astropy.coordinates import Angle, SkyCoord
from astropy.nddata import Cutout2D
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

def rms(array):
	sq_arr = np.square(array)
	avg = np.nanmean(sq_arr)
	return np.sqrt(avg)

def mask(reg, cutout):#masks everything except the region
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int') 

def plot_grid(datacube, masks, reject, snr_list, idx):
    n_images = len(datacube)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots + 1
    gs1 = gs.GridSpec(yplots, xplots, wspace=0.0, hspace=0.0, top=1.-0.5/(xplots+1), bottom=0.5/(xplots+1), left=0.5/(yplots+1), right=1-0.5/(yplots+1))
    plt.figure(figsize=(10, 0))
    for i in range(n_images):
        if reject[i] == True:
            colmap = 'gray'
        else:
            colmap = 'viridis'
        image = datacube[i]
        plt.subplot(gs1[i])
        plt.imshow(image, origin='lower', cmap=colmap)
        plt.imshow(masks[i], origin='lower', cmap='gray', alpha=0.2)
        plt.xticks([])
        plt.yticks([])
        plt.text(0.22,0.2,'id '+str(idx[i])+' snr '+str(np.around(snr_list[i],2)),{'size':7,'color':'white'})

def reject_sources(name, catname, datname, min_snr=5, max_size=None, max_size_ID=None, flux_type = 'peak'):

	#catname: name of catalog of sources, with required columns 'gauss_x_'+name, 'gauss_y_'+name, 'FWHM_major_'+name, 'FWHM_minor_'+name, and 'position_angle_'+name
	#datname: name of data fits image
	#min_snr: minimum SNR value, all sources with SNR below this will be rejected
	#max_size: maximum major axis radius for ellipse around source, in sigma
	#max_size_ID: alternatively, give the index of a source to set that source's radius to the maximum
	#flux_type: if 'peak', chooses the brightest pixel, if 'percentile', flux measured by average of top 10% of pixels

	catalog = fits.getdata(catname)
	catalog = Table(catalog)

	bad_inds = np.where(np.isnan(catalog['ap_flux_'+name])==True)
	catalog.remove_rows(bad_inds)

	fl = fits.open(datname)
	data = fl[0].data.squeeze()
	header = fl[0].header
	mywcs = WCS(header).celestial

	sigma_to_FWHM = 2*np.sqrt(2*np.log(2))

	pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg #for conversion to pixels

	snr_vals = []
	cutout_images = []
	masks = []
	bg_arr = []
	bg_arr2 = []
	reject = np.full(len(catalog), False)

	for i in range(len(catalog)):
		x_cen = catalog['gauss_x_'+name][i] * u.deg
		y_cen = catalog['gauss_y_'+name][i] * u.deg
		major_fwhm = (catalog['FWHM_major_'+name][i] * u.arcsec).to(u.degree)
		minor_fwhm = (catalog['FWHM_minor_'+name][i] * u.arcsec).to(u.degree)
		position_angle = catalog['position_angle_'+name][i] * u.deg

		annulus_width = 15
		center_pad = 10 #pad between ellipse and inner radius
		# Define some ellipse properties in pixel coordinates
		position = SkyCoord(x_cen, y_cen, frame='icrs', unit=(u.deg,u.deg))
		pix_position = np.array(position.to_pixel(mywcs))

		pix_major_fwhm = major_fwhm/pixel_scale
		pix_minor_fwhm = minor_fwhm/pixel_scale

		# Cutout section of the image we care about, to speed up computation time
		size = ((center_pad+annulus_width)*pixel_scale+major_fwhm)*2.2 #2.2 is arbitrary to get entire annulus and a little extra
		cutout = Cutout2D(data, position, size, mywcs, mode='partial') #cutout of outer circle
		cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1]) #center of the cutout in pixel coords

		# Define the aperture regions needed for SNR
		ellipse_reg = regions.EllipsePixelRegion(cutout_center, pix_major_fwhm*2., pix_minor_fwhm*2., angle=position_angle)
		innerann_reg = regions.CirclePixelRegion(cutout_center, center_pad+pix_major_fwhm)
		outerann_reg = regions.CirclePixelRegion(cutout_center, center_pad+pix_major_fwhm+annulus_width)

		# Make masks from aperture regions
		ellipse_mask = mask(ellipse_reg, cutout)
		annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

		# Calculate the SNR and aperture flux sums
		pixels_in_annulus = cutout.data[annulus_mask.astype('bool')] #pixels within annulus
		pixels_in_ellipse = cutout.data[ellipse_mask.astype('bool')] #pixels in ellipse
		bg_rms = rms(pixels_in_annulus)
		bg_mean = np.mean(pixels_in_annulus)
		bg_median = np.median(pixels_in_annulus)

		if flux_type == 'peak':
			peak_flux = catalog['peak_flux_'+name][i]
		if flux_type == 'percentile':
			top_percent = np.nanpercentile(pixels_in_ellipse, 90)
			peak_flux = np.mean(pixels_in_ellipse[pixels_in_ellipse > top_percent])
		snr = peak_flux / bg_rms
		catalog['ap_flux_err_'+name][i] = bg_rms
		bg_arr.append(bg_mean)
		bg_arr2.append(bg_median)
		

		if snr < min_snr: #if low snr, reject
			reject[i] = True
		if max_size_ID is not None:
			if catalog['major_sigma'][i] >  catalog['major_sigma'][max_size_ID] + 0.01/3600: #if too big a source, reject
				reject[i] = True
		if max_size is not None:
			if catalog['major_sigma'][i] >  max_size: #if too big a source, reject
				reject[i] = True
		snr_vals.append(snr)	
		cutout_images.append(cutout.data)
		masks.append(ellipse_mask + annulus_mask)

	catalog['bg_mean_'+name] = bg_arr
	catalog['bg_median_'+name] = bg_arr2
	plot_grid(cutout_images, masks, reject, snr_vals, catalog['_idx_'+name])
	plt.show(block=False)
	line_remove = input('enter id values for sources to exclude from the catalog, seperated by whitespace: ')
	man_input_rem = np.array(line_remove.split(), dtype='int')
	id_ind_true = np.where(np.in1d(catalog['_idx_'+name],man_input_rem) == True)
	reject[id_ind_true] = True

	line_keep = input('enter id values for removed sources to include from the catalog, seperated by whitespace: ')
	man_input_keep = np.array(line_keep.split(), dtype='int')
	id_ind_false = np.where(np.in1d(catalog['_idx_'+name],man_input_keep) == True)
	reject[id_ind_false] = False

	rej_ind = np.where(reject == True)
	catalog.remove_rows(rej_ind)
	return catalog
