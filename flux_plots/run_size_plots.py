import matplotlib.pyplot as plt
import numpy as np

from size_plots import *

from astropy.io import fits
from astropy.table import Table
from functools import reduce

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/r0.5_catalog_best.fits'))

B3_deconv_maj = data['fwhm_maj_deconv_B3'][np.where(data['fit_flag_B3']=='y')[0]]
B3_deconv_maj_err = data['fwhm_maj_deconv_err_B3'][np.where(data['fit_flag_B3']=='y')[0]]

B6_deconv_maj = data['fwhm_maj_deconv_B6'][np.where(data['fit_flag_B6']=='y')[0]]
B6_deconv_maj_err = data['fwhm_maj_deconv_err_B6'][np.where(data['fit_flag_B6']=='y')[0]]

B7_deconv_maj = data['fwhm_maj_deconv_B7'][np.where(data['fit_flag_B7']=='y')[0]]
B7_deconv_maj_err = data['fwhm_maj_deconv_err_B7'][np.where(data['fit_flag_B7']=='y')[0]]

#disk_size_hist([B3_deconv_maj, B6_deconv_maj, B7_deconv_maj], ['B3', 'B6', 'B7'])

disk_size_hist([B3_deconv_maj, B6_deconv_maj, B7_deconv_maj], ['B3 fwhm major', 'B6 fwhm major', 'B7 fwhm major'], 'size_hist_allbands.png')

B3B6ind = np.intersect1d(np.where(data['fit_flag_B3']=='y')[0], np.where(data['fit_flag_B6']=='y')[0])
B3_deconv_maj_B3B6ind = data['fwhm_maj_deconv_B3'][B3B6ind]
B6_deconv_maj_B3B6ind = data['fwhm_maj_deconv_B6'][B3B6ind]

disk_size_hist([B3_deconv_maj_B3B6ind - B6_deconv_maj_B3B6ind], ['B3 - B6 fwhm major deconv'], 'B3B6_difference_size_hist.png')
