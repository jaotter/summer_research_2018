from astropy.table import Table
from lifelines import KaplanMeierFitter
import astropy.constants as constants
from astropy.modeling import blackbody
from KM_plot import plot_KM
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def calc_dmass(fluxes, freq, dist):
    Tdust = 20*u.K
    kappa0 = 2*u.cm**2/u.g
    nu0 = constants.c/(1.3*u.mm)

    Bnu = blackbody.blackbody_nu(freq, 20*u.K)
    Dmass = (fluxes*dist**2)/(kappa0*(freq/nu0)*Bnu)
    Dmass = (Dmass.decompose()).to(u.earthMass*u.sr)

    return Dmass

tab_path = '/home/jotter/nrao/tables'

dmass_data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_may21_calc_vals_mask.fits')

eis_data = Table.read(f'{tab_path}/eisner_tbl.txt', format='ascii')
lupus_data = Table.read(f'{tab_path}/LupusDisks_Ansdell2016_dist_combined.txt', format='ascii', data_start=1)
ophi_data = Table.read(f'{tab_path}/Ophiuchus_Williams2019.txt', format='ascii')
perseus_data = Table.read(f'{tab_path}/Perseus_Anderson2019.txt', format='ascii', header_start=2, data_start=4, data_end=63, delimiter='\t')
taurus_data = Table.read(f'{tab_path}/TaurusDisks_Andrews2005.txt', format='ascii')
sco_data = Table.read(f'{tab_path}/UpperSco_Barenfield2016.txt', format='ascii')

IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')

nonIR_src = np.setdiff1d(dmass_data['ID'], IR_tab['ID'])
nonIR_ind = [np.where(dmass_data['ID']==d_id)[0][0] for d_id in nonIR_src]
IR_ind = [np.where(dmass_data['ID']==d_id)[0][0] for d_id in IR_tab['ID']]
omc1 = dmass_data[nonIR_ind]
onc = dmass_data[IR_ind]

IR_nondet = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_nondet_may21_full_ulim.fits')
B3_ulims = IR_nondet['B3_flux_ulim'].data * u.mJy

#all in units of earthMass
B3_mdust1 = dmass_data['dust_mass_B3'].data
#B3_ulims = np.repeat(0.05, 52)*u.mJy
#B3_ulims = np.repeat(0.05, 230)*u.mJy
B3_mdust_ulim = calc_dmass(B3_ulims, 98*u.GHz, 400*u.pc).value
B3_mdust = np.concatenate((B3_mdust1, B3_mdust_ulim))
B3_mdust_flag = np.concatenate((np.repeat(True, len(B3_mdust1)), np.repeat(False, len(B3_mdust_ulim))))

omc1_B3mdust = omc1['dust_mass_B3'].data
omc1_mdust_flag = np.repeat(True, len(omc1_B3mdust))

onc_B3mdust = onc['dust_mass_B3'].data
onc_mdust_flag = np.repeat(True, len(onc_B3mdust))

B3_mdust_onc = np.concatenate((onc_B3mdust, B3_mdust_ulim))
B3_mdust_flag_onc = np.concatenate((np.repeat(True, len(onc_B3mdust)), np.repeat(False, len(B3_mdust_ulim))))

#B3_mdust_onc = onc_B3mdust
#B3_mdust_flag_onc = np.repeat(True, len(onc_B3mdust))

lupus_mdust = lupus_data['MDust'].data
lupus_mdust_flag = np.repeat(True, len(lupus_mdust))

sco_mdust = sco_data['Mdust'].data
sco_mdust_flag = np.where(sco_data['f_Mdust'] == '<', False, True)

eis_mdust_str = eis_data['M_dust^a'].data
eis_mdust1 = np.array([float(mdust.split(' ')[0]) for mdust in eis_mdust_str])
eis_mdust_flag1 = np.repeat(True, len(eis_mdust1))
#eis_mdust = eis_mdust[eis_mdust>0]
eis_nondet = Table.read('/home/jotter/nrao/tables/eisner_nondetect.txt', format='ascii', delimiter='\t')#, data_start=2, header_start=1, 
print(eis_nondet.info)
eis_ulim = eis_nondet['F_lambda850mum']
eis_ulim_flux = []
for ul in eis_ulim:
    eis_ulim_flux.append(float(ul[1:]))
eis_freq = (constants.c/(850*u.micron)).to(u.GHz)
eis_mdust2 = calc_dmass(eis_ulim_flux*u.mJy, eis_freq, 400*u.pc).value
eis_mdust_flag2 = np.repeat(False, len(eis_mdust2))
eis_mdust = np.concatenate((eis_mdust1, eis_mdust2))
eis_mdust_flag = np.concatenate((eis_mdust_flag1, eis_mdust_flag2))



ophiucus_flux = ophi_data['F225'].data * u.mJy
ophi_freq = 225*u.GHz
dists = ophi_data['d'].data*u.pc
ophi_mdust = calc_dmass(ophiucus_flux, ophi_freq, dists).value
#ophi_mdust2 = 0.58*ophiucus_flux #conversion from paper, uses same equation
#ophi_mdust = ophi_mdust[ophi_mdust>0]

ophi_mdust_flag = np.repeat(True, len(ophi_mdust))

taurus_flux = taurus_data['F450'].data*u.mJy
lamb = 450*u.micron
freq = constants.c/lamb
taurus_mdust = calc_dmass(taurus_flux, freq, 140*u.pc).value
taurus_mdust_flag = np.repeat(True, len(taurus_mdust))
taurus_mdust_flag[np.where(taurus_data['l_F450'] == '<')] = False
#taurus_mdust_flag = taurus_mdust_flag[taurus_mdust>0]
#taurus_mdust = taurus_mdust[taurus_mdust>0]

#print(B3_mdust)

#B3_mdust = B3_mdust[B3_mdust < 300]


#plot_KM([eis_mdust, lupus_mdust, sco_mdust, B3_mdust, B6_mdust, B7_mdust, ophi_mdust, taurus_mdust], ['E18', 'Lupus', 'Upper Sco', 'B3', 'B6', 'B7', 'Ophiucus', 'Taurus'],
#        [eis_mdust_flag, lupus_mdust_flag, sco_mdust_flag, B3_mdust_flag, B6_mdust_flag, B7_mdust_flag, ophi_mdust_flag, taurus_mdust_flag],
#        savepath='/home/jotter/nrao/plots/KM_dust_mass.pdf')

#plot_KM([eis_mdust, lupus_mdust, sco_mdust, ophi_mdust, B3_mdust_onc, taurus_mdust], ['E18', 'Lupus', 'Upper Sco', 'Ophiucus', 'ONC B3', 'Taurus'],
#        [eis_mdust_flag, lupus_mdust_flag, sco_mdust_flag, ophi_mdust_flag, B3_mdust_flag_onc, taurus_mdust_flag], savepath='/home/jotter/nrao/plots/KM_dust_mass_aug20_onc_noulim.pdf', left_censor=True, cdf=False)

onc_combined = np.concatenate((B3_mdust_onc, eis_mdust))
onc_combined_flag = np.concatenate((B3_mdust_flag_onc, eis_mdust_flag))


plot_KM([lupus_mdust, sco_mdust, ophi_mdust, B3_mdust_onc, taurus_mdust], ['Lupus', 'Upper Sco', 'Ophiucus', 'ONC (B3)', 'Taurus'],
        [lupus_mdust_flag, sco_mdust_flag, ophi_mdust_flag, B3_mdust_flag_onc, taurus_mdust_flag], savepath='/home/jotter/nrao/plots/KM_dust_mass_may21_onc.pdf', left_censor=True, cdf=False, plot_quantity='Mdust')
