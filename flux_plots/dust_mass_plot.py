from astropy.table import Table
from lifelines import KaplanMeierFitter

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, savepath='/home/jotter/nrao/summer_research_2018/plots/KM_comparisons.png'):

    kmf = KaplanMeierFitter()

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    for ind in range(len(arrays)):
        print(labels[ind])
        kmf.fit(arrays[ind], label=labels[ind])
        kmf.plot(ax=ax)

    plt.legend()
    ax.set_xlabel(r'$\log(M_{dust}/M_{\oplus})$')
    ax.set_ylabel(r'$P \geq M_{dust}$')
    plt.savefig(savepath)


tab_path = '/home/jotter/nrao/tables'

eis_data = Table.read(f'{tab_path}/eisner_tbl.txt', format='ascii')
dmass_data = Table.read(f'{tab_path}/inf_vals_all_updt.fits')
lupus_data = Table.read(f'{tab_path}/LupusDisks_Ansdell2016.txt', format='ascii')
ophi_data = Table.read(f'{tab_path}/Ophiucus_FluxSize_Cieza2018.txt', format='ascii.fixed_width')
perseus_data = Table.read(f'{tab_path}/Perseus_Anderson2019.txt', format='ascii', header_start=2, data_start=4, data_end=63, delimiter='\t')
taurus_data = Table.read(f'{tab_path}/TaurusDisks_Andrews2005.txt', format='ascii')
sco_data = Table.read(f'{tab_path}/UpperSco_Barenfield2016.txt', format='ascii')

#all in units of earthMass
eis_mdust = eis_data['M_dust^a'].data
lupus_mdust = lupus_data['MDust'].data
sco_mdust = sco_data['Mdust'].data
B3_mdust = dmass_data['dust_mass_B3'].data
B6_mdust = dmass_data['dust_mass_B6'].data
B7_mdust = dmass_data['dust_mass_B7'].data

print(sco_mdust)




#check Ophiucus data - need to convert flux to disk mass, same with taurus

perseus_mdisk = perseus_data['M_disk^c'] #also total disk mass I think
taurus_mdisk = taurus_data['Mass'] #this is total disk mass - not comparable


plot_KM([eis_mdust, lupus_mdust, sco_mdust, B3_mdust], ['E18', 'Lupus', 'Upper Sco', 'B3'])
