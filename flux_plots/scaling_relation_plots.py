from astropy.table import Table
from astropy.io import fits, ascii
from radio_beam import Beam

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

import astropy.constants as constants

band='B3'

#three plots: R and L, R and M, L and M
data = Table.read('../tables/table_meas_B3.fits')
data_inf = Table.read('../tables/inf_vals_all.fits')

andrews = ascii.read('../tables/andrews_table.txt', format='tab')
andrews.rename_column('$\mathrm{log}{R}_{\mathrm{eff}}/\mathrm{au}$', 'R_eff')
andrews.rename_column('F_nu/mJy', 'Fnu')
andrews.rename_column('$\mathrm{log}{L}_{\mathrm{mm}}/\mathrm{Jy}$', 'L_mm')
andrews.remove_row(0)

andrews_table4 = ascii.read('../tables/andrews_table4.txt', format='tab')
andrews_table4.rename_column('$d/\mathrm{pc}$', 'd')

R_eff_andrw = []
R_eff_err_up_andrw = []
R_eff_err_down_andrw = []

Fnu_andrw = []
Lmm_scaled_andrw = []

dists_andrw = []

for ind in range(len(andrews)):
    if len(andrews['R_eff'][ind]) > 10:
        R_eff_andrw.append(float(andrews['R_eff'][ind][0:4]))
        R_eff_err_up_andrw.append(float(andrews['R_eff'][ind][-6:-3]))
        R_eff_err_down_andrw.append(float(andrews['R_eff'][ind][10:14]))
        
        Fnu_andrw.append(float(andrews['Fnu'][ind].split('$')[0]))
        Lmm_scaled_andrw.append(float(andrews['L_mm'][ind][0:4]))
        dists_andrw.append(float(andrews_table4['d'][ind][0:3]))

Lmm_scaled_andrw = 10**np.array(Lmm_scaled_andrw)*u.Jy
R_eff_andrw = 10**np.array(R_eff_andrw)*u.au
Fnu_andrw = np.array(Fnu_andrw)*u.mJy
dists_andrw = np.array(dists_andrw)*u.pc

andrews_freq = 335*u.GHz
andrews_beam = Beam(major=0.29*u.arcsec, minor=0.25*u.arcsec, pa=89*u.degree)
R_andrw = (np.sin(0.25*u.arcsec)*dists_andrw).to(u.au)/2

#T_B = L_mm.to(u.K, equivalencies=u.brightness_temperature(andrews_freq))
#val = ((2*constants.h*(andrews_freq**3))/(L_mm*(constants.c**2))).decompose()
#T_B = ((constants.h*andrews_freq/constants.k_B)*(np.log(1+val))**(-1)).decompose()
T_B_andrw = Fnu_andrw.to(u.K, andrews_beam.jtok_equiv(andrews_freq))
L_mm_andrw = (4*np.pi*(R_andrw**2)*constants.sigma_sb*(T_B_andrw)**4).to(u.L_sun)

fig = plt.figure(figsize=(5,5))

deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_'+band]) == False)[0]

Rarcsec = data['fwhm_maj_deconv_'+band][deconv_ind]*u.arcsec
d = (414*u.pc).to(u.au)
Rau = (Rarcsec.to(u.rad)*d).value
Rau_err = (((data['fwhm_maj_err_'+band][deconv_ind])*u.arcsec).to(u.rad)*d).value

scaled_B3flux = data['ap_flux_B3'][deconv_ind]*(414/140)**2

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

coeff = ((((2-q)*x*constants.c**2)/(4*np.pi*(andrews_freq**2)*constants.k_B*scriptF*T0*(r0**q)))**(1/(2-q)))
#coeff = (coeff.value*u.m).to(u.au)

Lstar1 = 1
Lstar10 = 10
Lmm_arr = np.logspace(-3, 0, 4)*u.Jy

R_eqn9_1 = (Lstar1**(-1/(4*(2-q)))*Lmm_arr**(1/(2-q))*coeff/SIGMA_TO_FWHM).decompose()
R_eqn9_10 = (Lstar10**(-1/(4*(2-q)))*Lmm_arr**(1/(2-q))*coeff/SIGMA_TO_FWHM).decompose()



srcIind = np.where(data['D_ID'][deconv_ind] == 10)[0]
BNind = np.where(data['D_ID'][deconv_ind] == 20)[0]
print('source I scaled Lmm: %f' % (scaled_B3flux[srcIind]))
print('source BN scaled Lmm: %f' % (scaled_B3flux[BNind]))

print('source I R: %f' % (Rau[srcIind]))
print('source BN R: %f' % (Rau[BNind]))


plt.errorbar(scaled_B3flux, Rau, yerr=Rau_err, linestyle='', marker='o', label='Band 3 measurements')
plt.errorbar(Lmm_scaled_andrw.value, FWHM_andrw, yerr=FWHM_err_andrw, linestyle='', marker='o', label='Andrews et al. 2018')

#plt.plot(Lmm_arr, R_eqn9_1, linestyle='-', marker='', label='A18 Eqn 9, '+r'$L_* = 1 L_\odot$')
#plt.plot(Lmm_arr, R_eqn9_10, linestyle='-', marker='', label='A18 Eqn 9, '+r'$L_* = 10 L_\odot$')
plt.legend()
plt.xlabel('Scaled luminosity (Jy)')
plt.ylabel('R (au)')
plt.yscale('log')
plt.xscale('log')

plt.savefig('plots/scaling_rels_scaledflux.png', dpi=300)

'''
plt.errorbar(Lmm, Rau, yerr=Rau_err, linestyle='', marker='o', label='Band 3 measurements')
plt.errorbar(L_mm_andrw.value, FWHM_andrw, yerr=FWHM_err_andrw, linestyle='', marker='o', label='Andrews et al. 2018')
plt.legend()
plt.xlabel('L lower limit (Lsun)')
plt.ylabel('R (au)')
plt.yscale('log')
plt.xscale('log')

plt.savefig('plots/scaling_rels_flux.png')


plt.clf()
plt.figure()
flux_B3_hist, bins = np.histogram(data['ap_flux_B3'], density=False, bins=20)
andrews_flux_hist, bins2  = np.histogram(Fnu_andrw.value/1000, bins=bins, density=False)

print(data['ap_flux_B3'])
print(Fnu_andrw.value/1000)

plotpts = []
widths = []
for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))
                                             

plt.bar(plotpts, flux_B3_hist, widths, label='B3 data', alpha=0.5)
plt.bar(plotpts, andrews_flux_hist, widths, label='Andrews et. al. 2018', alpha=0.5)
plt.legend()
plt.savefig('plots/andrews_flux_hist.png',dpi=300)
'''
