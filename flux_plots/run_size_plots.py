import matplotlib.pyplot as plt
import numpy as np

from size_plots import *

from astropy.io import fits
from astropy.table import Table
from functools import reduce

data = Table(fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits'))

B3_deconv_maj = data['fwhm_maj_deconv_B3']
B3_deconv_maj_err = data['fwhm_maj_deconv_err_B3']
B6_deconv_maj = data['fwhm_maj_deconv_B6']
B6_deconv_maj_err = data['fwhm_maj_deconv_err_B6']
B7_deconv_maj = data['fwhm_maj_deconv_B7']
B7_deconv_maj_err = data['fwhm_maj_deconv_err_B7']

B3deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0]

#R_hist_eisner(data['fwhm_maj_deconv_B3'][B3deconv_ind], 'B3', 'B3 size', 'size_hist_eisner.png')

size_comp_eisner('eisner_size_comp.png')

#disk_size_hist_3panel([data['fwhm_maj_B3'], data['fwhm_maj_B6'], data['fwhm_maj_B7']], ['B3', 'B6', 'B7'], 'size_hist_3panel_meas.png')

#disk_size_hist_3panel([B3_deconv_maj, B6_deconv_maj, B7_deconv_maj], ['B3', 'B6', 'B7'], 'size_hist_3panel_deconv.png')

ind13 = np.where(data['D_ID']==13)[0]

B6B7ind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_B7']) == False)[0], np.where(np.isnan(data['fwhm_maj_B6'])==False)[0])

B6B7ind = np.delete(B6B7ind, np.where(B6B7ind == ind13))

B6B7deconvind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_deconv_B7']) == False)[0], np.where(np.isnan(data['fwhm_maj_deconv_B6'])==False)[0])

B6B7deconvind = np.delete(B6B7deconvind, np.where(B6B7deconvind == ind13))
B6B7deconvind = np.delete(B6B7deconvind, np.where(B6B7deconvind == 3))

B6B7ind = np.setdiff1d(B6B7ind, B6B7deconvind)

#size_comp([data['fwhm_maj_B6'][B6B7ind], data['fwhm_maj_B7'][B6B7ind]], [data['fwhm_maj_deconv_B6'][B6B7deconvind], data['fwhm_maj_deconv_B7'][B6B7deconvind]],[data['fwhm_maj_err_B6'][B6B7ind], data['fwhm_maj_err_B7'][B6B7ind]], [data['fwhm_maj_deconv_err_B6'][B6B7deconvind], data['fwhm_maj_deconv_err_B7'][B6B7deconvind]], ['B6 fwhm', 'B7 fwhm'], 'size_comp_B6B7_maj.png')

B3B6ind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_B6'])==False)[0])

B3B6ind = np.delete(B3B6ind, np.where(B3B6ind == ind13))

B3B6deconvind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_deconv_B6'])==False)[0])

B3B6ind = np.setdiff1d(B3B6ind, B3B6deconvind)

#size_comp([data['fwhm_maj_B3'][B3B6ind], data['fwhm_maj_B6'][B3B6ind]], [data['fwhm_maj_deconv_B3'][B3B6deconvind], data['fwhm_maj_deconv_B6'][B3B6deconvind]],[data['fwhm_maj_err_B3'][B3B6ind], data['fwhm_maj_err_B6'][B3B6ind]], [data['fwhm_maj_deconv_err_B6'][B3B6deconvind], data['fwhm_maj_deconv_err_B6'][B3B6deconvind]], ['B3 fwhm', 'B6 fwhm'], 'size_comp_B3B6_maj.png')

B3B7ind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_B7'])==False)[0])

B3B7ind = np.delete(B3B7ind, np.where(B3B7ind == ind13))

B3B7deconvind = np.intersect1d(np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0], np.where(np.isnan(data['fwhm_maj_deconv_B7'])==False)[0])

B3B6deconvind = np.delete(B3B6deconvind, np.where(B3B6deconvind == ind13))
B3B6deconvind = np.delete(B3B6deconvind, np.where(B3B6deconvind == 3))

B3B7ind = np.setdiff1d(B3B7ind, B3B7deconvind)

#size_comp_simple([data['fwhm_maj_deconv_B3'][B3B6deconvind], data['fwhm_maj_deconv_B6'][B3B6deconvind]], [data['fwhm_maj_deconv_err_B3'][B3B6deconvind], data['fwhm_maj_deconv_err_B6'][B3B6deconvind]], ['B3 deconvolved FWHM (as)', 'B6 deconvolved FWHM (as)'], 'size_comp_simple_B3B6.png')

#size_comp_simple([data['fwhm_maj_deconv_B6'][B6B7deconvind], data['fwhm_maj_deconv_B7'][B6B7deconvind]], [data['fwhm_maj_deconv_err_B6'][B6B7deconvind], data['fwhm_maj_deconv_err_B7'][B6B7deconvind]], ['B6 deconvolved FWHM (as)', 'B7 deconvolved FWHM (as)'], 'size_comp_simple_B6B7.png')

#size_comp_simple([data['fwhm_maj_deconv_B3'][B3B7deconvind], data['fwhm_maj_deconv_B7'][B3B7deconvind]], [data['fwhm_maj_deconv_err_B3'][B3B7deconvind], data['fwhm_maj_deconv_err_B7'][B3B7deconvind]], ['B3 deconvolved FWHM (as)', 'B7 deconvolved FWHM (as)'], 'size_comp_simple_B3B7.png')
