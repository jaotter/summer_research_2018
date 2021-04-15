from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Circle
from mpl_plot_templates import asinh_norm
from matplotlib.lines import Line2D
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/Forbrich2016_r0.5_nov20.fits')

fig = plt.figure(figsize=(5,5))

B3flux = np.log10(tab['ap_flux_B3']*1000)
B3flux_nonlog = tab['ap_flux_B3']


plt.plot(B3flux, np.log10(tab['Spk']), marker='*', linestyle='')
plt.ylabel(r'$\log(S_{5GHz})$ (mJy)')
plt.xlabel(r'$\log(S_{98GHz})$ (mJy)')

plt.savefig('/home/jotter/nrao/plots/vla_flux_plot.png', bbox_inches='tight')

plt.clf()


calc = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_nov20_calc_vals.fits')
calc_ind = []
for ind in range(len(tab)):
    calc_ind.append(np.where(calc['D_ID'] == tab[ind]['D_ID'])[0])

calc = calc[np.array(calc_ind)]
alpha_calc = (calc['alpha_B3B6']+calc['alpha_B6B7'])/2
nanind = np.where(np.isnan(alpha_calc) == False)[0]

vla_alpha = tab['alpha'][nanind]
calc_alpha = alpha_calc[nanind]

fig2 = plt.figure(figsize=(5,5))

plt.plot(calc_alpha, vla_alpha, marker='o', linestyle='')
plt.ylabel(r'$\alpha_{5GHz}$')
plt.xlabel(r'$\alpha_{\sim 200GHz}$')
plt.xlim(-0.4, 1.75)
plt.ylim(-0.4, 1.75)
plt.savefig('/home/jotter/nrao/plots/vla_alpha_plot.png', bbox_inches='tight')

plt.clf()

pred_flux = np.log10(tab['Spk'])*((98/5)**(-0.1))
pred_flux_nonlog = tab['Spk']*((98/5)**(-0.1))


flux_diff = B3flux - pred_flux


fig = plt.figure(figsize=(5,5))

plt.plot(B3flux, flux_diff, marker='*', linestyle='')
plt.ylabel(r'$\log(S_{98GHz}) - \log(S_{predicted, 98GHz})$ (mJy)')
plt.xlabel(r'$\log(S_{98GHz})$ (mJy)')

plt.savefig('/home/jotter/nrao/plots/vla_flux_pred_plot_diff.png', bbox_inches='tight')

plt.clf()

fig = plt.figure(figsize=(5,5))

plt.plot(B3flux, pred_flux, marker='*', linestyle='')
plt.ylabel(r'$\log(S_{predicted, 98GHz})$ (mJy)')
plt.xlabel(r'$\log(S_{98GHz})$ (mJy)')
plt.xlim(-1.5,2)
plt.ylim(-1.5,2)
plt.savefig('/home/jotter/nrao/plots/vla_flux_pred_plot.png', bbox_inches='tight')

plt.clf()

##spatial plot

cube_path = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
fl = fits.open(cube_path)
img = fl[0].data.squeeze()
header = fl[0].header
mywcs = WCS(header).celestial

fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.15,0.1,0.8,0.8],projection=mywcs)
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
dec.set_major_formatter('dd:mm:ss.s')
ra.ticklabels.set_fontsize(18)
dec.ticklabels.set_fontsize(18)

ax.imshow(img, origin='lower', transform=ax.get_transform(mywcs), norm=asinh_norm.AsinhNorm(), vmin=-0.001, vmax=0.005, cmap='Greys')


flux_ratio = B3flux_nonlog/pred_flux_nonlog

ff_ind = np.where(flux_ratio <= 0.5)[0]
med_ind = np.where(flux_ratio < 1.5)[0]
int_ind = np.setdiff1d(med_ind, ff_ind)
dust_ind = np.where(flux_ratio >= 1.5)[0]

print(len(flux_ratio), len(ff_ind), len(int_ind), len(dust_ind))

B3_pix = mywcs.all_world2pix(tab['RA_B3']*u.degree, tab['DEC_B3']*u.degree, 0)

for ind in ff_ind:
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=100, fill=False, color='tab:red')
    ax.add_patch(circ)
    #ax.text(B3_pix[0][ind]-1, B3_pix[1][ind]+3, B3_table['D_ID'][ind], color='red')

for ind in int_ind:
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=100, fill=False, color='tab:green')
    ax.add_patch(circ)
    #ax.text(B3_pix[0][ind]-1, B3_pix[1][ind]+3, B3_table['D_ID'][ind], color='red')

for ind in dust_ind:
    circ = Circle((B3_pix[0][ind], B3_pix[1][ind]), radius=100, fill=False, color='tab:blue')
    ax.add_patch(circ)
    #ax.text(B3_pix[0][ind]-1, B3_pix[1][ind]+3, B3_table['D_ID'][ind], color='red')

legend_elements = [Line2D([0], [0], color='tab:red', marker='o', fillstyle='none', linestyle=''),
                   Line2D([0], [0], color='tab:green', marker='o', fillstyle='none', linestyle=''),
                   Line2D([0], [0], color='tab:blue', marker='o', fillstyle='none', linestyle='')]

ax.legend(legend_elements, [r'$S_{98GHz}/S_{ff} \leq 0.5$', r'$0.5 < S_{98GHz}/S_{ff} < 1.5$', r'$S_{98GHz}/S_{ff} \geq 1.5$'])

    
plt.savefig('/home/jotter/nrao/plots/vla_overlay_plot_nonlog.png', bbox_inches='tight')
