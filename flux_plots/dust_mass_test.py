from astropy.modeling import blackbody
from astropy.table import Table
import astropy.units as u
import numpy as np
import astropy.constants as constants
import matplotlib.pyplot as plt

eis = Table.read('../tables/eisner_tbl.txt', format='ascii')
eis_Fdust = []
for fd in eis['F_dust']:
    eis_Fdust.append(float(fd.split()[0]))

data = Table.read('../tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')
B3flux = data['ap_flux_B3']
B3freq = 98*u.GHz

B7flux = data['ap_flux_B7']

for i in B7flux:
    print(i*1000)
print(eis_Fdust)

Tdust = 20*u.K
kappa0 = 2*u.cm**2/u.g
dist = 414*u.pc
nu0 = (constants.c/(1.3*u.mm)).decompose()
nu = (constants.c/(850*u.um)).decompose()
Bnu = blackbody.blackbody_nu(nu, 20*u.K)
#Bnu = Bnu.decompose().to(u.J/u.m**2/u.sr)

B3nu = 98*u.GHz
Bnu_B3 = blackbody.blackbody_nu(B3nu, 20*u.K)

B3_Dmass = ((B3flux*u.Jy)*dist**2)/(kappa0*(B3nu/nu0)*Bnu_B3)
B3_Dmass = np.log10(B3_Dmass.decompose().to(u.earthMass*u.sr).value)

eis_Dmass = ((eis_Fdust*u.mJy)*dist**2)/(kappa0*(nu/nu0)*Bnu)
eis_Dmass = np.log10(eis_Dmass.decompose().to(u.earthMass*u.sr).value)
eis_Dmass = eis_Dmass[np.where(np.isinf(eis_Dmass) == False)[0]]
    
all_hist, bins_B3 = np.histogram(np.concatenate((B3_Dmass, eis_Dmass)))
B3hist, b = np.histogram(B3_Dmass, bins=bins_B3)
eishist, b = np.histogram(eis_Dmass, bins=bins_B3)

plotpts = []
widths = []

for b in range(len(bins_B3[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(bins_B3[b] + (bins_B3[b+1]-bins_B3[b])/2)
    widths.append((bins_B3[b+1]-bins_B3[b]))
                    
plt.bar(plotpts, B3hist, widths, label='B3', alpha=0.5)
plt.bar(plotpts, eishist, widths, label='eisner', alpha=0.5)
plt.legend()
plt.savefig('plots/Dmass_hist2.png')

inf_vals = Table.read('../tables/inf_vals_all.fits')

