##script to calculate where sources are significantly different in diff bands


from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u

import matplotlib.pyplot as plt
import numpy as np

data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')

maj_B3 = data['fwhm_maj_deconv_B3']
maj_B3_err = data['fwhm_maj_deconv_err_B3']

maj_B7 = data['fwhm_maj_deconv_B7']
maj_B7_err = data['fwhm_maj_deconv_err_B7']

maj_B6 = data['fwhm_maj_deconv_B6']
maj_B6_err = data['fwhm_maj_deconv_err_B6']


print('B3 vs B6')
bigger_ind_B3B6 = np.where(maj_B3 - 3*maj_B3_err > maj_B6 + 3*maj_B6_err)[0]
print('Bigger in B3 than B6 ', data['D_ID'][bigger_ind_B3B6])

smaller_ind_B3B6 = np.where(maj_B3 + 3*maj_B3_err < maj_B6 - 3*maj_B6_err)[0]
print('Bigger in B6 than B3 ', data['D_ID'][smaller_ind_B3B6])


print('B3 vs B7')
bigger_ind_B3B7 = np.where(maj_B3 - 3*maj_B3_err > maj_B7 + 3*maj_B7_err)[0]
print('Bigger in B3 than B7 ', data['D_ID'][bigger_ind_B3B7])

smaller_ind_B3B7 = np.where(maj_B7 - 3*maj_B7_err > maj_B3 + 3*maj_B3_err)[0]
print('Bigger in B7 than B3 ', data['D_ID'][smaller_ind_B3B7])


print('B6 vs B7')
bigger_ind_B6B7 = np.where(maj_B6 - 3*maj_B6_err > maj_B7 + 3*maj_B7_err)[0]
print('Bigger in B6 than B7 ', data['D_ID'][bigger_ind_B6B7])

smaller_ind_B6B7 = np.where(maj_B7 - 3*maj_B7_err > maj_B6 + 3*maj_B6_err)[0]
print('Bigger in B7 than B6 ', data['D_ID'][smaller_ind_B6B7])
