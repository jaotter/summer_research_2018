from astropy.table import Table
from lifelines import KaplanMeierFitter

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, upper_lim_flags, savepath='/home/jotter/nrao/plots/KM_dust_mass_onlyB3.pdf'):

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
    ax.set_xlabel(r'$\log(M_{dust}/M_{\oplus})$')
    ax.set_ylabel(r'$P \leq M_{dust}$')
    ax.set_xlim(np.log10(0.01), 4)
    plt.savefig(savepath)


tab_path = '/home/jotter/nrao/tables'

dmass_data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_jun20_calc_vals.fits')

eis_data = Table.read(f'{tab_path}/eisner_tbl.txt', format='ascii')
lupus_data = Table.read(f'{tab_path}/LupusDisks_Ansdell2016.txt', format='ascii')
ophi_data = Table.read(f'{tab_path}/Ophiucus_FluxSize_Cieza2018.txt', format='ascii.fixed_width', delimiter=' ', data_start=2)
perseus_data = Table.read(f'{tab_path}/Perseus_Anderson2019.txt', format='ascii', header_start=2, data_start=4, data_end=63, delimiter='\t')
taurus_data = Table.read(f'{tab_path}/TaurusDisks_Andrews2005.txt', format='ascii')
sco_data = Table.read(f'{tab_path}/UpperSco_Barenfield2016.txt', format='ascii')

#all in units of earthMass
B3_mdust = np.log10(dmass_data['dust_mass_B3'].data)
B6_mdust = np.log10(dmass_data['dust_mass_B6'].data)
B7_mdust = np.log10(dmass_data['dust_mass_B7'].data)
B3_mdust_flag = np.repeat(True, len(B3_mdust))

lupus_mdust = np.log10(lupus_data['MDust'].data)
lupus_mdust_flag = np.repeat(True, len(lupus_mdust))

sco_mdust = np.log10(sco_data['Mdust'].data)
sco_mdust_flag = np.where(sco_data['f_Mdust'] == '<', False, True)

eis_mdust_str = eis_data['M_dust^a'].data
eis_mdust = np.array([float(mdust.split(' ')[0]) for mdust in eis_mdust_str])
eis_mdust = np.log10(eis_mdust[eis_mdust > 0])
eis_mdust_flag = np.repeat(True, len(eis_mdust))

ophiucus_flux = ophi_data['F1.3'].data
ophi_mdust = 0.58*ophiucus_flux #conversion from paper, uses same equation
ophi_mdust_flag = np.repeat(True, len(ophi_mdust))

perseus_mdisk = perseus_data['M_disk^c'] #also total disk mass I think

taurus_mdust = taurus_data['Mass'].data*u.Msun
taurus_mdust = np.log10(taurus_mdust.to(u.earthMass).value)
taurus_mdust_flag = np.repeat(True, len(taurus_mdust))
#print(B3_mdust)

#B3_mdust = B3_mdust[B3_mdust < 300]

B6_mdust = B6_mdust[np.isnan(B6_mdust)==False]
B6_mdust_flag = np.repeat(True, len(B6_mdust))
B7_mdust = B7_mdust[np.isnan(B7_mdust)==False]
B7_mdust_flag = np.repeat(True, len(B7_mdust))

#plot_KM([eis_mdust, lupus_mdust, sco_mdust, B3_mdust, B6_mdust, B7_mdust, ophi_mdust, taurus_mdust], ['E18', 'Lupus', 'Upper Sco', 'B3', 'B6', 'B7', 'Ophiucus', 'Taurus'],
#        [eis_mdust_flag, lupus_mdust_flag, sco_mdust_flag, B3_mdust_flag, B6_mdust_flag, B7_mdust_flag, ophi_mdust_flag, taurus_mdust_flag],
#        savepath='/home/jotter/nrao/plots/KM_dust_mass.pdf')
plot_KM([eis_mdust, lupus_mdust, sco_mdust, B3_mdust, ophi_mdust, taurus_mdust], ['E18', 'Lupus', 'Upper Sco', 'B3', 'Ophiucus', 'Taurus'],
        [eis_mdust_flag, lupus_mdust_flag, sco_mdust_flag, B3_mdust_flag, ophi_mdust_flag, taurus_mdust_flag])
