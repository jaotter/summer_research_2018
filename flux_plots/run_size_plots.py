import matplotlib.pyplot as plt
import numpy as np

from size_plots import *

from astropy.io import fits
from astropy.table import Table
from functools import reduce

data = Table(fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted.fits'))

B3_deconv_maj = data['fwhm_maj_deconv_B3']
B3_deconv_maj_err = data['fwhm_maj_deconv_err_B3']
B6_deconv_maj = data['fwhm_maj_deconv_B6']
B6_deconv_maj_err = data['fwhm_maj_deconv_err_B6']
B7_deconv_maj = data['fwhm_maj_deconv_B7']
B7_deconv_maj_err = data['fwhm_maj_deconv_err_B7']

B6B7ind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_B7']) == False)[0], np.where(np.isnan(data['fwhm_maj_B6'])==False)[0])

B6B7ind = np.delete(B6B7ind, np.where(B6B7ind == 12))

B6B7deconvind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_deconv_B7']) == False)[0], np.where(np.isnan(data['fwhm_maj_deconv_B6'])==False)[0])

B6B7ind = np.setdiff1d(B6B7ind, B6B7deconvind)

#size_comp([data['fwhm_maj_B6'][B6B7ind], data['fwhm_maj_B7'][B6B7ind]], [data['fwhm_maj_deconv_B6'][B6B7deconvind], data['fwhm_maj_deconv_B7'][B6B7deconvind]],[data['fwhm_maj_err_B6'][B6B7ind], data['fwhm_maj_err_B7'][B6B7ind]], [data['fwhm_maj_deconv_err_B6'][B6B7deconvind], data['fwhm_maj_deconv_err_B7'][B6B7deconvind]], ['B6 fwhm', 'B7 fwhm'], 'size_comp_B6B7_maj.png')

B3B6ind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_B6'])==False)[0])

B3B6ind = np.delete(B3B6ind, np.where(B3B6ind == 12))

B3B6deconvind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_deconv_B6'])==False)[0])

B3B6ind = np.setdiff1d(B3B6ind, B3B6deconvind)

size_comp([data['fwhm_maj_B3'][B3B6ind], data['fwhm_maj_B6'][B3B6ind]], [data['fwhm_maj_deconv_B3'][B3B6deconvind], data['fwhm_maj_deconv_B6'][B3B6deconvind]],[data['fwhm_maj_err_B3'][B3B6ind], data['fwhm_maj_err_B6'][B3B6ind]], [data['fwhm_maj_deconv_err_B6'][B3B6deconvind], data['fwhm_maj_deconv_err_B6'][B3B6deconvind]], ['B3 fwhm', 'B6 fwhm'], 'size_comp_B3B6_maj.png')

B3B7ind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_B7'])==False)[0])

B3B7ind = np.delete(B3B7ind, np.where(B3B7ind == 12))

B3B7deconvind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_deconv_B7'])==False)[0])

B3B7ind = np.setdiff1d(B3B7ind, B3B7deconvind)

size_comp([data['fwhm_maj_B3'][B3B7ind], data['fwhm_maj_B7'][B3B7ind]], [data['fwhm_maj_deconv_B3'][B3B7deconvind], data['fwhm_maj_deconv_B7'][B3B7deconvind]],[data['fwhm_maj_err_B3'][B3B7ind], data['fwhm_maj_err_B7'][B3B7ind]], [data['fwhm_maj_deconv_err_B7'][B3B7deconvind], data['fwhm_maj_deconv_err_B7'][B3B7deconvind]], ['B3 fwhm', 'B7 fwhm'], 'size_comp_B3B7_maj.png')


#disk_size_hist([B3_deconv_maj, B6_deconv_maj, B7_deconv_maj], ['B3 fwhm major', 'B6 fwhm major', 'B7 fwhm major'], 'size_hist_allbands_deconv.png')

#disk_size_hist([data['fwhm_maj_B3'], data['fwhm_maj_B6'], data['fwhm_maj_B7']], ['B3 fwhm major', 'B6 fwhm major', 'B7 fwhm major'], 'size_hist_allbands.png')
