from astropy.table import Table
from astropy.io import fits, ascii

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

band='B3'

#three plots: R and L, R and M, L and M
data = Table.read('../tables/measured_vals_all.fits')
data_inf = Table.read('../tables/inf_vals_all.fits')

andrews = ascii.read('../tables/andrews_table.txt', format='tab')
andrews.rename_column('$\mathrm{log}{R}_{\mathrm{eff}}/\mathrm{au}$', 'R_eff')
andrews.rename_column('$\mathrm{log}{L}_{\mathrm{mm}}/\mathrm{Jy}$', 'L_mm')
andrews.remove_row(0)

R_eff = []
R_eff_err = []

L_mm = []
L_mm_err = []

for ind in range(len(andrews)):
    if len(andrews['R_eff'][ind]) > 10:
        R_eff.append(float(andrews['R_eff'][ind][0:4]))
        R_err_up = float(andrews['R_eff'][ind][-6:-3])
        R_err_down = float(andrews['R_eff'][ind][10:14])
        R_eff_err.append((R_err_up, R_err_down))

        L_mm.append(float(andrews['L_mm'][ind][0:4]))
        L_err_up = float(andrews['L_mm'][ind][-6:-3])
        L_err_down = float(andrews['L_mm'][ind][10:14])
        L_mm_err.append((L_err_up, L_err_down))

print(R_eff, L_mm)
R_eff = 10**R_eff
L_mm = 10**L_mm

fig = plt.figure(figsize=(5,5))

deconv_ind = np.where(np.isnan(data['fwhm_maj_deconv_'+band]) == False)[0]

Rarcsec = data['fwhm_maj_deconv_'+band][deconv_ind]*u.arcsec
d = (414*u.pc).to(u.au)
Rau = (Rarcsec.to(u.rad)*d).value
Rau_err = (((data['fwhm_maj_err_'+band][deconv_ind])*u.arcsec).to(u.rad)*d).value

Lmm = data_inf['lower_lum_'+band][deconv_ind]

Mdust = data_inf['dust_mass_'+band][deconv_ind]
Mdust_err = data_inf['dust_mass_err_'+band][deconv_ind]


plt.errorbar(Lmm, Rau, yerr=Rau_err, linestyle='', marker='o')
plt.xlabel('L lower limit (Lsun)')
plt.ylabel('R (au)')
plt.yscale('log')
plt.xscale('log')

plt.savefig('plots/scaling_rels.png')
