import matplotlib.pyplot as plt
import numpy as np
import sys
from size_plots import *

from astropy.io import fits
from astropy.table import Table
from functools import reduce

data = Table(fits.getdata('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_feb21_ulim.fits'))

B3_deconv_maj = data['fwhm_maj_deconv_B3']
B3_deconv_maj_err = data['fwhm_maj_deconv_err_B3']
B6_deconv_maj = data['fwhm_maj_deconv_B6']
B6_deconv_maj_err = data['fwhm_maj_deconv_err_B6']
B7_deconv_maj = data['fwhm_maj_deconv_B7']
B7_deconv_maj_err = data['fwhm_maj_deconv_err_B7']

B3deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0]

#R_hist_eisner(data['fwhm_maj_deconv_B3'][B3deconv_ind], 'B3 size', 'size_hist_eisner.pdf')

size_comp_eisner('eisner_size_comp.pdf')

#disk_size_hist_3panel([data['fwhm_maj_deconv_B3'], data['fwhm_maj_deconv_B6'], data['fwhm_maj_deconv_B7']], ['B3', 'B6', 'B7'], 'size_hist_3panel_ulim.pdf', ulim_arrs=[data['upper_lim_B3'], data['upper_lim_B6'], data['upper_lim_B7']])


IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_feb21_full.fits')
nonIR_src = np.setdiff1d(data['Seq'], IR_tab['Seq'])
nonIR_ind = [np.where(data['Seq']==d_id)[0][0] for d_id in nonIR_src]
IR_src = IR_tab['Seq']
IR_ind = [np.where(data['Seq']==d_id)[0][0] for d_id in IR_src]

#disk_size_hist_3panel_IR([B3_deconv_maj[IR_ind], B6_deconv_maj[IR_ind], B7_deconv_maj[IR_ind]], [B3_deconv_maj[nonIR_ind], B6_deconv_maj[nonIR_ind], B7_deconv_maj[nonIR_ind]], ['B3', 'B6', 'B7'], 'size_hist_3panel_deconv_IR.pdf')

B3_deconv_ONC = B3_deconv_maj[IR_ind]
B3_deconv_ind_ONC = np.where(np.isnan(B3_deconv_ONC) == False)[0]

B3_deconv_OMC1 = B3_deconv_maj[nonIR_ind]
B3_deconv_ind_OMC1 = np.where(np.isnan(B3_deconv_OMC1) == False)[0]


#R_hist_eisner(B3_deconv_ONC[B3_deconv_ind_ONC], 'B3 ONC sources', 'size_hist_eisner_all.pdf', norm=False, size_arr2=B3_deconv_OMC1[B3_deconv_ind_OMC1], label2='B3 OMC1 sources')
#R_hist_eisner(B3_deconv_ONC[B3_deconv_ind_ONC], 'B3 ONC sources', 'size_hist_eisner_ONC.pdf', norm=False)
#R_hist_eisner(B3_deconv_OMC1[B3_deconv_ind_OMC1], 'B3 OMC1 sources', 'size_hist_eisner_OMC1.pdf', norm=False)

disk_size_hist_3panel([B3_deconv_maj, B6_deconv_maj, B7_deconv_maj], ['B3', 'B6', 'B7'], 'size_hist_3panel_deconv.pdf')

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

#size_comp_simple([data['fwhm_maj_deconv_B3'][B3B6deconvind], data['fwhm_maj_deconv_B6'][B3B6deconvind]], [data['fwhm_maj_deconv_err_B3'][B3B6deconvind], data['fwhm_maj_deconv_err_B6'][B3B6deconvind]], ['B3 deconvolved FWHM (as)', 'B6 deconvolved FWHM (as)'], 'size_comp_simple_B3B6.pdf', src_ids = data['D_ID'][B3B6deconvind])

#size_comp_simple([data['fwhm_maj_deconv_B6'][B6B7deconvind], data['fwhm_maj_deconv_B7'][B6B7deconvind]], [data['fwhm_maj_deconv_err_B6'][B6B7deconvind], data['fwhm_maj_deconv_err_B7'][B6B7deconvind]], ['B6 deconvolved FWHM (as)', 'B7 deconvolved FWHM (as)'], 'size_comp_simple_B6B7.pdf', src_ids = data['D_ID'][B6B7deconvind])

#size_comp_simple([data['fwhm_maj_deconv_B3'][B3B7deconvind], data['fwhm_maj_deconv_B7'][B3B7deconvind]], [data['fwhm_maj_deconv_err_B3'][B3B7deconvind], data['fwhm_maj_deconv_err_B7'][B3B7deconvind]], ['B3 deconvolved FWHM (as)', 'B7 deconvolved FWHM (as)'], 'size_comp_simple_B3B7.pdf', src_ids = data['D_ID'][B3B7deconvind])
