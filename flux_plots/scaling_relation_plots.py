from astropy.table import Table
from astropy.io import fits, ascii
from radio_beam import Beam

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

import astropy.constants as constants

from scipy.stats import linregress

band='B3'

#three plots: R and L, R and M, L and M
data = Table.read('../tables/r0.5_catalog_bgfit_apr20.fits')
data_inf = Table.read('../tables/r0.5_apr20_calc_vals.fits')

andrews = ascii.read('/home/jotter/nrao/tables/andrews_table.txt', format='tab')
andrews.rename_column('$\mathrm{log}{R}_{\mathrm{eff}}/\mathrm{au}$', 'R_eff')
andrews.rename_column('F_nu/mJy', 'Fnu')
andrews.rename_column('$\mathrm{log}{L}_{\mathrm{mm}}/\mathrm{Jy}$', 'L_mm')
andrews.remove_row(0)

andrews_table4 = ascii.read('/home/jotter/nrao/tables/andrews_table4.txt', format='tab')
andrews_table4.rename_column('$d/\mathrm{pc}$', 'd')

R_eff_andrw_log = []
R_eff_err_up_andrw = []
R_eff_err_down_andrw = []

Fnu_andrw = []
Lmm_scaled_andrw_log = []
Lmm_scaled_err_up_andrw = []
Lmm_scaled_err_down_andrw = []
dists_andrw = []

for ind in range(len(andrews)):
    if len(andrews['R_eff'][ind]) > 10:
        R_eff_andrw_log.append(float(andrews['R_eff'][ind][0:4]))
        R_eff_err_up_andrw.append(float(andrews['R_eff'][ind][-6:-2]))
        R_eff_err_down_andrw.append(float(andrews['R_eff'][ind][9:14]))
        
        Fnu_andrw.append(float(andrews['Fnu'][ind].split('$')[0]))
        Lmm_scaled_andrw_log.append(float(andrews['L_mm'][ind][0:5]))
        Lmm_scaled_err_up_andrw.append(float(andrews['L_mm'][ind][-6:-2]))
        Lmm_scaled_err_down_andrw.append(float(andrews['L_mm'][ind][10:15]))
        
        dists_andrw.append(float(andrews_table4['d'][ind][0:3]))


Lmm_scaled_andrw = 10**np.array(Lmm_scaled_andrw_log)*u.Jy
Lmm_scaled_err_up_andrw = 10**np.array(np.array(Lmm_scaled_andrw_log)+np.array(Lmm_scaled_err_up_andrw))*u.Jy - Lmm_scaled_andrw
Lmm_scaled_err_down_andrw = 10**np.array(np.array(Lmm_scaled_andrw_log)-np.array(Lmm_scaled_err_down_andrw))*u.Jy - Lmm_scaled_andrw
Lmm_scaled_err_andrw = np.array((Lmm_scaled_err_up_andrw.value, Lmm_scaled_err_down_andrw.value))
R_eff_andrw = 10**np.array(R_eff_andrw_log)*u.au
Fnu_andrw = np.array(Fnu_andrw)*u.mJy
dists_andrw = np.array(dists_andrw)*u.pc

R_eff_err_up_andrw = 10**np.array(np.array(R_eff_andrw_log) + np.array(R_eff_err_up_andrw))*u.au - R_eff_andrw
R_eff_err_down_andrw = 10**np.array(np.array(R_eff_andrw_log) - np.array(R_eff_err_down_andrw))*u.au - R_eff_andrw

andrews_freq = 335*u.GHz

'''andrews_beam = Beam(major=0.29*u.arcsec, minor=0.25*u.arcsec, pa=89*u.degree)
R_andrw = (np.sin(0.25*u.arcsec)*dists_andrw).to(u.au)/2

#T_B = L_mm.to(u.K, equivalencies=u.brightness_temperature(andrews_freq))
#val = ((2*constants.h*(andrews_freq**3))/(L_mm*(constants.c**2))).decompose()
#T_B = ((constants.h*andrews_freq/constants.k_B)*(np.log(1+val))**(-1)).decompose()
T_B_andrw = Fnu_andrw.to(u.K, andrews_beam.jtok_equiv(andrews_freq))
L_mm_andrw = (4*np.pi*(R_andrw**2)*constants.sigma_sb*(T_B_andrw)**4).to(u.L_sun)'''

fig = plt.figure(figsize=(5,5))

deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_'+band]) == False)[0]

Rarcsec = data['fwhm_maj_deconv_'+band][deconv_ind]*u.arcsec
d = (414*u.pc).to(u.au)
Rau = (Rarcsec.to(u.rad)*d).value
Rau_err = (((data['fwhm_maj_err_'+band][deconv_ind])*u.arcsec).to(u.rad)*d).value

scaled_B3flux = data['ap_flux_B3'][deconv_ind]*(414/140)**2
scaled_B3flux_err = data['ap_flux_err_B3'][deconv_ind]*(414/140)**2

Lmm = data_inf['lower_lum_'+band][deconv_ind]

Mdust = data_inf['dust_mass_'+band][deconv_ind]
Mdust_err = data_inf['dust_mass_err_'+band][deconv_ind]


R_eff_err_andrw = np.array((R_eff_err_up_andrw, R_eff_err_down_andrw))

SIGMA_TO_FWHM = np.sqrt(8*np.log(2))

FWHM_andrw = R_eff_andrw.value*SIGMA_TO_FWHM
FWHM_err_andrw = R_eff_err_andrw*SIGMA_TO_FWHM

plt.figure()


q = 0.57
x = 0.68
scriptF = 0.3
T0 = 30*u.K
r0 = 10*u.au

Rau_arr = np.logspace(0,3,5)*u.au
Lmm_eqn9_01 = ((Rau_arr*SIGMA_TO_FWHM)**(2-q)*andrews_freq**2*(.1)**(.25)*constants.k_B*scriptF*T0*r0**q/((2-q)*constants.c**2*(140*u.pc)**2)).decompose().to(u.Jy)
Lmm_eqn9 = ((Rau_arr*SIGMA_TO_FWHM)**(2-q)*andrews_freq**2*constants.k_B*scriptF*T0*r0**q/((2-q)*constants.c**2*(140*u.pc)**2)).decompose().to(u.Jy)
Lmm_eqn9_10 = ((Rau_arr*SIGMA_TO_FWHM)**(2-q)*andrews_freq**2*(10)**(.25)*constants.k_B*scriptF*T0*r0**q/((2-q)*constants.c**2*(140*u.pc)**2)).decompose().to(u.Jy)


srcIind = np.where(data['D_ID'][deconv_ind] == 30)[0]
BNind = np.where(data['D_ID'][deconv_ind] == 43)[0]
#print('source I scaled Lmm: %f' % (scaled_B3flux[srcIind]))
#print('source BN scaled Lmm: %f' % (scaled_B3flux[BNind]))

#print('source I R: %f' % (Rau[srcIind]))
#print('source BN R: %f' % (Rau[BNind]))

scaled_B3flux_rem = np.delete(scaled_B3flux, [srcIind[0], BNind[0]])
Rau_rem = np.delete(Rau, [srcIind[0], BNind[0]])
Rau_rem_err = np.delete(Rau_err, [srcIind[0], BNind[0]])

linreg_params_B3 = linregress(np.log10(scaled_B3flux_rem), np.log10(Rau_rem))
Lmm_arr = np.logspace(-3, 0, 5)
linreg_line_B3 = np.log10(Lmm_arr)*linreg_params_B3[0] + linreg_params_B3[1]

#plt.plot(Lmm_arr, 10**linreg_line_B3, linestyle='-', marker='', label='Fit to B3 data')
print('linreg params for B3 fit', linreg_params_B3)

linreg_params_andrw = linregress(np.log10(Lmm_scaled_andrw.value), np.log10(FWHM_andrw))
linreg_line_andrw = np.log10(Lmm_arr)*linreg_params_andrw[0] + linreg_params_andrw[1]

#plt.plot(Lmm_arr, 10**linreg_line_andrw, linestyle='-', marker='', label='Fit to A18 data')
print('linreg params for A18 fit', linreg_params_andrw)

linreg_params_all = linregress(np.log10(np.concatenate((Lmm_scaled_andrw.value, scaled_B3flux_rem))), np.log10(np.concatenate((FWHM_andrw, Rau_rem))))
linreg_line_all = np.log10(Lmm_arr)*linreg_params_all[0] + linreg_params_all[1]

print(len(Lmm_scaled_andrw.value) + len(scaled_B3flux_rem))

plt.plot(Lmm_arr, 10**linreg_line_all, linestyle='-', marker='', label='Fit to B3 and A18 data')
print('linreg params for all fit', linreg_params_all)

plt.errorbar(scaled_B3flux, Rau, yerr=Rau_err, xerr=scaled_B3flux_err, linestyle='', marker='o', label='Band 3 measurements')
plt.errorbar(Lmm_scaled_andrw.value, FWHM_andrw, yerr=FWHM_err_andrw, xerr=Lmm_scaled_err_andrw, linestyle='', marker='o', label='Andrews et al. 2018')

#plt.plot(Lmm_eqn9_01, Rau_arr, linestyle='-', marker='', label='A18 Eqn 9, '+r'$L_* = 0.1 L_\odot$')
plt.plot(Lmm_eqn9, Rau_arr, linestyle='--', marker='', label='A18 Eqn 9, '+r'$L_* = 1 L_\odot$')
#plt.plot(Lmm_eqn9_10, Rau_arr, linestyle='-', marker='', label='A18 Eqn 9, '+r'$L_* = 10 L_\odot$')
plt.legend()
plt.xlim(0.001, 1)
plt.ylim(7, 600)
plt.xlabel('Scaled luminosity (Jy)')
plt.ylabel('R (AU)')
plt.yscale('log')
plt.xscale('log')

#plt.savefig('/home/jotter/nrao/plots/scaling_rels_scaledflux_allfit.png', dpi=300)

fig = plt.figure()

d = 414*u.pc
d = d.to(u.AU)

size_arr = data['fwhm_maj_deconv_B3']

size_arr = size_arr[np.isnan(size_arr)==False]*u.arcsec
size_arr = (size_arr.to(u.rad)*d).value
size_arr = np.log10(size_arr)

eisner_data = ascii.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='tab')
eisner_ind = np.where(eisner_data['R_disk'] != '<5')[0]
eisner_R = [float(x.split()[0])*2 for x in eisner_data['R_disk'][eisner_ind]]
eisner_R = np.log10(eisner_R)

FWHM_andrw = np.log10(FWHM_andrw)

total_R = np.concatenate((FWHM_andrw, eisner_R, size_arr))
total_hist_R, bins_R = np.histogram(total_R, bins=10)
andrw_hist, bins = np.histogram(FWHM_andrw, density=False, bins=bins_R)
eis_hist, bins = np.histogram(eisner_R, density=False, bins=bins_R)
hist, bins = np.histogram(size_arr, density=False, bins=bins_R)


plotpts = []
widths = []
for b in range(len(bins[:-1])): #creating points to plot - midpoints of bbins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))


plt.bar(plotpts, hist, widths, edgecolor = 'black', label='band 3 sizes', alpha=0.4)
plt.bar(plotpts, eis_hist, widths, edgecolor='black', label='E18 sizes', alpha=0.4)
plt.bar(plotpts, andrw_hist, widths, edgecolor='black', label='A18 sizes', alpha=0.4)

plt.xlabel('log(Major FWHM / AU)')
plt.ylabel('Number')
plt.legend()

plt.savefig('/home/jotter/nrao/plots/R_hist_all.png', dpi=500)


fig = plt.figure()
scaled_B3flux = np.log10(scaled_B3flux)

Lmm_scaled_andrw = np.log10(Lmm_scaled_andrw.value)

L_tot = np.concatenate((scaled_B3flux, Lmm_scaled_andrw))
hist_Ltot, bins_L = np.histogram(L_tot, density=False, bins=10)

hist_L, bins = np.histogram(scaled_B3flux, density=False, bins=bins_L)
andrw_hist_L, bins = np.histogram(Lmm_scaled_andrw, density=False, bins=bins_L)

plotpts = []
widths = []

for b in range(len(bins[:-1])): #creating points to plot - midpoints of bbins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))

plt.bar(plotpts, hist_L, widths, edgecolor = 'black', label='band 3 sizes', alpha=0.4)
plt.bar(plotpts, andrw_hist_L, widths, edgecolor='black', label='A18 sizes', alpha=0.4)

plt.xlabel('log(Scaled Luminosity / Jy)')
plt.ylabel('Number')
plt.legend()

plt.savefig('/home/jotter/nrao/plots/L_hist_all.png', dpi=500)
