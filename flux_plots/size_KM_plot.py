from astropy.table import Table
from lifelines import KaplanMeierFitter

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, upper_lim_flags, savepath='/home/jotter/nrao/plots/KM_size_plot.pdf'):

    kmf = KaplanMeierFitter()

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    for ind in range(len(arrays)):
        print(labels[ind])
        if upper_lim_flags[ind] is not None:
            kmf.fit_left_censoring(arrays[ind], upper_lim_flags[ind], label=labels[ind])
        else:
            kmf.fit(arrays[ind], upper_lim_flags[ind], label=labels[ind])
        
        kmf.plot(ax=ax)

    plt.legend()
    ax.set_xlabel(r'FWHM (AU)')
    ax.set_ylabel(r'$P \leq$ FWHM')
    #ax.set_xlim(np.log10(0.01), 4)
    plt.savefig(savepath)


tab_path = '/home/jotter/nrao/tables'

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_jun20.fits')

eis_data = Table.read(f'{tab_path}/eisner_tbl.txt', format='ascii')
lupus_data = Table.read(f'{tab_path}/LupusDisks_Ansdell2016_dist_combined.txt', format='ascii')
ophi_data = Table.read(f'{tab_path}/Ophiucus_FluxSize_Cieza2018.txt', format='ascii.fixed_width', delimiter=' ', data_start=2)
#perseus_data = Table.read(f'{tab_path}/Perseus_Anderson2019.txt', format='ascii', header_start=2, data_start=4, data_end=63, delimiter='\t')
#taurus_data = Table.read(f'{tab_path}/TaurusDisks_Andrews2005.txt', format='ascii')
sco_data = Table.read(f'{tab_path}/UpperSco_Barenfeld2016_size.txt', format='ascii')

#all in units of AU


B3_fwhm = data['fwhm_maj_deconv_B3']*u.arcsec
B3_fwhm = B3_fwhm[np.where(np.isnan(B3_fwhm)==False)]
B3_fwhm_au = ((B3_fwhm.to(u.radian))*(415*u.pc).to(u.AU)).value

B6_fwhm = data['fwhm_maj_deconv_B6'].data*u.arcsec
B5_fwhm = B6_fwhm[np.where(np.isnan(B6_fwhm)==False)]
B6_fwhm_au = ((B6_fwhm.to(u.radian))*(415*u.pc).to(u.AU)).value

B7_fwhm = data['fwhm_maj_deconv_B7'].data*u.arcsec
B7_fwhm = B7_fwhm[np.where(np.isnan(B7_fwhm)==False)]
B7_fwhm_au = ((B7_fwhm.to(u.radian))*(415*u.pc).to(u.AU)).value

B3_fwhm_flag = np.repeat(True, len(B3_fwhm))

lupus_fwhm = lupus_data['a']*u.arcsec
lupus_fwhm = lupus_fwhm[np.where(lupus_fwhm != 0)]
lupus_dist = lupus_data['Dis']*u.pc
lupus_dist = lupus_dist[np.where(lupus_fwhm != 0)]
lupus_fwhm_au = (lupus_fwhm.to(u.radian)*lupus_dist.to(u.AU)).value
lupus_fwhm_flag = np.repeat(True, len(lupus_fwhm_au))

sco_fwhm = sco_data['FWHM'][np.where(np.isnan(sco_data['FWHM'])==False)]
sco_fwhm_au = (sco_fwhm.to(u.radian)*(145*u.pc).to(u.AU)).value
sco_fwhm_flag = np.repeat(True, len(sco_fwhm)) #np.where(sco_data['f_Mdust'] == '<', False, True)

eis_rdisk_str = eis_data['R_disk'].data
eis_ulim_ind = np.where(eis_rdisk_str == '<5')
eis_rdisk_str[eis_ulim_ind] = '5 '
eis_fwhm = np.array([float(rdisk.split(' ')[0]) for rdisk in eis_rdisk_str])*2
eis_fwhm_flag = np.repeat(True, len(eis_fwhm))
eis_fwhm_flag[eis_ulim_ind] = False

ulim_ophi_ind = np.where(ophi_data['Major'] == '...')
ophi_fwhm = ophi_data['Major']
ophi_fwhm[ulim_ophi_ind] = 200
ophi_fwhm = np.array(ophi_fwhm,dtype='float')*u.mas
ophi_fwhm_au = (ophi_fwhm.to(u.radian)*(140*u.pc).to(u.AU)).value
ophi_fwhm_flag = np.repeat(True, len(ophi_fwhm_au))
ophi_fwhm_flag[ulim_ophi_ind] = False

plot_KM([eis_fwhm, lupus_fwhm_au, sco_fwhm_au, B3_fwhm_au, ophi_fwhm_au], ['E18', 'Lupus', 'Upper Sco', 'B3', 'Ophiucus'],
        [eis_fwhm_flag, lupus_fwhm_flag, sco_fwhm_flag, B3_fwhm_flag, ophi_fwhm_flag])
