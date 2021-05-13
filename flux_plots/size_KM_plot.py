from astropy.table import Table
from lifelines import KaplanMeierFitter

from KM_plot import plot_KM, KM_median, bootstrap_ci
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
tab_path = '/home/jotter/nrao/tables'

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')

eis_data = Table.read(f'{tab_path}/eisner_tbl.txt', format='ascii')
lupus_data = Table.read(f'{tab_path}/LupusDisks_Ansdell2016_dist_combined.txt', format='ascii')
ophi_data = Table.read(f'{tab_path}/Ophiucus_FluxSize_Cieza2018.txt', format='ascii.fixed_width', delimiter=' ', data_start=2)
#perseus_data = Table.read(f'{tab_path}/Perseus_Anderson2019.txt', format='ascii', header_start=2, data_start=4, data_end=63, delimiter='\t')
#taurus_data = Table.read(f'{tab_path}/TaurusDisks_Andrews2005.txt', format='ascii')

#all in units of AU


B3_fwhm = data['fwhm_maj_deconv_B3']*u.arcsec
B3_fwhm_flag = np.repeat(True, len(B3_fwhm))
B3_ulim_ind = np.where(np.isnan(data['upper_lim_B3'])==False)
B3_fwhm_flag[B3_ulim_ind] = False
B3_fwhm = ((B3_fwhm.to(u.radian))*(400*u.pc).to(u.AU)).value
B3_fwhm[B3_ulim_ind] = data['upper_lim_B3'][B3_ulim_ind].data
B3_hwhm_au = B3_fwhm/2

B3_deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B3']) == False)[0]
B3_fwhm_resolved = B3_fwhm[B3_deconv_ind]
B3_fwhm_flag_resolved = np.repeat(True, len(B3_fwhm_resolved))

print('B3')
median, lower, upper = KM_median(B3_fwhm, B3_fwhm_flag, return_type='percentile')
#ci = bootstrap_ci(100, B3_fwhm, B3_fwhm_flag)
print(f'{median} - {lower} + {upper}')

#print('B3 resolved')
#median = KM_median(B3_fwhm_resolved, B3_fwhm_flag_resolved, return_type='percentile', left_censor=False)
#ci = bootstrap_ci(100, B3_fwhm_resolved, B3_fwhm_flag_resolved)
#print(f'{median} +/- {ci}')

B6ind = np.where(np.isnan(data['ap_flux_B6']) == False)
B6_fwhm = data['fwhm_maj_deconv_B6'][B6ind]*u.arcsec
B6_fwhm_flag = np.repeat(True, len(B6_fwhm))
B6_ulim_ind = np.where(np.isnan(data['upper_lim_B6'][B6ind])==False)
B6_fwhm_flag[B6_ulim_ind] = False
B6_fwhm = ((B6_fwhm.to(u.radian))*(400*u.pc).to(u.AU)).value
B6_fwhm[B6_ulim_ind] = data['upper_lim_B6'][B6ind][B6_ulim_ind].data
B6_hwhm_au = B6_fwhm/2

print('B6')
median, lower, upper = KM_median(B6_fwhm, B6_fwhm_flag, return_type='percentile')
print(f'{median} - {lower} + {upper}')

#median = KM_median(B6_fwhm, B6_fwhm_flag, return_type='median')
#ci = bootstrap_ci(100, B6_fwhm, B6_fwhm_flag)
#print(f'{median} +/- {ci}')

B6_deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B6'][B6ind]) == False)[0]
B6_fwhm_resolved = B6_fwhm[B6_deconv_ind]
B6_fwhm_flag_resolved = np.repeat(True, len(B6_fwhm_resolved))

print('B6 resolved')
#median = KM_median(B6_fwhm_resolved, B6_fwhm_flag_resolved, return_type='median', left_censor=False)
#ci = bootstrap_ci(100, B6_fwhm_resolved, B6_fwhm_flag_resolved)
#print(f'{median} +/- {ci}')

B7ind = np.where(np.isnan(data['ap_flux_B7']) == False)
B7_fwhm = data['fwhm_maj_deconv_B7'][B7ind]*u.arcsec
B7_fwhm_flag = np.repeat(True, len(B7_fwhm))
B7_ulim_ind = np.where(np.isnan(data['upper_lim_B7'][B7ind])==False)
B7_fwhm_flag[B7_ulim_ind] = False
B7_fwhm = ((B7_fwhm.to(u.radian))*(400*u.pc).to(u.AU)).value
B7_fwhm[B7_ulim_ind] = data['upper_lim_B7'][B7ind][B7_ulim_ind].data
B7_hwhm_au = B7_fwhm/2

print('B7')
median, lower, upper = KM_median(B7_fwhm, B7_fwhm_flag, return_type='percentile')
print(f'{median} - {lower} + {upper}')

#median = KM_median(B7_fwhm, B7_fwhm_flag, return_type='median')
#ci = bootstrap_ci(100, B7_fwhm, B7_fwhm_flag)
#print(f'{median} +/- {ci}')

B7_deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_B7'][B7ind]) == False)[0]
B7_fwhm_resolved = B7_fwhm[B7_deconv_ind]
B7_fwhm_flag_resolved = np.repeat(True, len(B7_fwhm_resolved))

#print('B7 resolved')
#median = KM_median(B7_fwhm_resolved, B7_fwhm_flag_resolved, return_type='median', left_censor=False)
#ci = bootstrap_ci(100, B7_fwhm_resolved, B7_fwhm_flag_resolved)
#print(f'{median} +/- {ci}')

IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
nonIR_src = np.setdiff1d(data['ID'], IR_tab['ID'])
nonIR_ind = [np.where(data['ID']==d_id)[0][0] for d_id in nonIR_src]
IR_ind = [np.where(data['ID']==d_id)[0][0] for d_id in IR_tab['ID']]

onc = B3_hwhm_au[IR_ind]
onc_flag = B3_fwhm_flag[IR_ind]
omc1 = B3_hwhm_au[nonIR_ind]
omc1_flag = B3_fwhm_flag[nonIR_ind]

print('ONC')
#median = KM_median(onc*2, onc_flag, return_type='median')
#ci = bootstrap_ci(10000, onc*2, onc_flag)
#print(f'{median} +/- {ci}')
median, lower, upper = KM_median(onc*2, onc_flag*2, return_type='percentile')
print(f'{median} - {lower} + {upper}')


print('OMC1')
#median = KM_median(omc1, omc1_flag, return_type='median')
#ci = bootstrap_ci(1000, omc1, omc1_flag)
median, lower, upper = KM_median(omc1*2, omc1_flag*2, return_type='percentile')
print(f'{median} - {lower} + {upper}')

#print(f'{median} +/- {ci}')

B6_IR = np.intersect1d(B6ind, IR_ind)
nonIR_src_B6 = np.setdiff1d(data['ID'][B6ind], IR_tab['ID'])
nonIR_ind_B6 = [np.where(data['ID'][B6ind]==d_id)[0][0] for d_id in nonIR_src_B6]
IR_tab_ind_B6 = [np.where(IR_tab['ID']==d_id)[0][0] for d_id in data['ID'][B6_IR]]
IR_ind_B6 = [np.where(data['ID'][B6ind]==d_id)[0][0] for d_id in IR_tab['ID'][IR_tab_ind_B6]]
omc1_B6 = B6_hwhm_au[nonIR_ind_B6]
onc_B6 = B6_hwhm_au[IR_ind_B6]
onc_flag_B6 = B6_fwhm_flag[IR_ind_B6]
omc1_flag_B6 = B6_fwhm_flag[nonIR_ind_B6]

B7_IR = np.intersect1d(B7ind, IR_ind)
nonIR_src_B7 = np.setdiff1d(data['ID'][B7ind], IR_tab['ID'])
nonIR_ind_B7 = [np.where(data['ID'][B7ind]==d_id)[0][0] for d_id in nonIR_src_B7]
IR_tab_ind_B7 = [np.where(IR_tab['ID']==d_id)[0][0] for d_id in data['ID'][B7_IR]]
IR_ind_B7 = [np.where(data['ID'][B7ind]==d_id)[0][0] for d_id in IR_tab['ID'][IR_tab_ind_B7]]
omc1_B7 = B7_hwhm_au[nonIR_ind_B7]
onc_B7 = B7_hwhm_au[IR_ind_B7]
onc_flag_B7 = B7_fwhm_flag[IR_ind_B7]
omc1_flag_B7 = B7_fwhm_flag[nonIR_ind_B7]


lupus_fwhm = lupus_data['a']*u.arcsec
lupus_fwhm = lupus_fwhm[np.where(lupus_fwhm != 0)]
lupus_dist = lupus_data['Dis']*u.pc
lupus_dist = lupus_dist[np.where(lupus_fwhm != 0)]
lupus_fwhm_au = (lupus_fwhm.to(u.radian)*lupus_dist.to(u.AU)).value
lupus_fwhm_flag = np.repeat(True, len(lupus_fwhm_au))

lupus_hwhm_au = lupus_fwhm_au/2

print('Lupus')
#median = KM_median(lupus_fwhm_au, lupus_fwhm_flag, return_type='median')
#ci = bootstrap_ci(10000, lupus_fwhm_au, lupus_fwhm_flag)
#print(f'{median} +/- {ci}')
median, lower, upper = KM_median(lupus_fwhm_au, lupus_fwhm_flag, return_type='percentile')
print(f'{median} - {lower} + {upper}')



sco_r = []
with open('/home/jotter/nrao/tables/UpperSco_Barenfeld2017_edit.txt', 'r') as fl:
    for line in fl:
        sco_r.append(float(line.split('\t')[2].split(' ')[0]))
#sco_fwhm = sco_data['FWHM'][np.where(np.isnan(sco_data['FWHM'])==False)]
#sco_fwhm_au = (sco_fwhm.to(u.radian)*(145*u.pc).to(u.AU)).value
sco_fwhm_flag = np.repeat(True, len(sco_r)) #np.where(sco_data['f_Mdust'] == '<', False, True)

eis_rdisk_str = eis_data['R_disk'].data
eis_ulim_ind = np.where(eis_rdisk_str == '<5')
eis_rdisk_str[eis_ulim_ind] = '5 '
eis_hwhm = np.array([float(rdisk.split(' ')[0]) for rdisk in eis_rdisk_str])
eis_hwhm_flag = np.repeat(True, len(eis_hwhm))
eis_hwhm_flag[eis_ulim_ind] = False

print('E18')
#median = KM_median(eis_hwhm*2, eis_hwhm_flag, return_type='median')
#ci = bootstrap_ci(1000, eis_hwhm*2, eis_hwhm_flag)
#print(f'{median} +/- {ci}')

median, lower, upper = KM_median(eis_hwhm*2, eis_hwhm_flag, return_type='percentile')
print(f'{median} - {lower} + {upper}')

print('ONC combined')
onc_combined = np.concatenate((onc, eis_hwhm))*2
onc_combined_flag = np.concatenate((onc_flag, eis_hwhm_flag))
#median = KM_median(onc_combined, onc_combined_flag, return_type='median')
#ci = bootstrap_ci(1000, onc_combined, onc_combined_flag)
#print(f'{median} +/- {ci}')
median, lower, upper = KM_median(onc_combined, onc_combined_flag, return_type='percentile')
print(f'{median} - {lower} + {upper}')


#eis_fwhm = np.delete(eis_fwhm, eis_ulim_ind)
#eis_fwhm_flag = np.repeat(True, len(eis_fwhm))

ulim_ophi_ind = np.where(ophi_data['Major'] == '...')
ophi_fwhm = ophi_data['Major']
ophi_fwhm[ulim_ophi_ind] = 200
ophi_fwhm = np.array(ophi_fwhm,dtype='float')*u.mas
ophi_fwhm_au = (ophi_fwhm.to(u.radian)*(140*u.pc).to(u.AU)).value
ophi_fwhm_flag = np.repeat(True, len(ophi_fwhm_au))
ophi_fwhm_flag[ulim_ophi_ind] = False

ophi_hwhm_au = ophi_fwhm_au/2

#ophi_fwhm_au = np.delete(ophi_fwhm_au, ulim_ophi_ind)
#ophi_fwhm_flag = np.repeat(True, len(ophi_fwhm_au))


#plot_KM([eis_hwhm, lupus_hwhm_au, sco_r, ophi_hwhm_au, onc, omc1], ['E18', 'Lupus', 'Upper Sco', 'Ophiuchus', 'ONC B3', 'OMC1 B3'],
#        [eis_hwhm_flag, lupus_fwhm_flag, sco_fwhm_flag, ophi_fwhm_flag, onc_flag, omc1_flag], savepath='/home/jotter/nrao/plots/KM_size_plot_may21_hwhm_omc1_onc.pdf', left_censor=True,
#        plot_quantity='Rdisk')
