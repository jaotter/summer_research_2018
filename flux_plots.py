from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits')
B3nu = (84+116)/2
B6nu = (211+275)/2

#B3/B6 plot:
ind = np.concatenate((np.where(np.isnan(data['flux_B6'])== False)[0],np.where(np.isnan(data['flux_B3'])==False)[0]))
B3_plot = data['flux_B3'][ind]
B6_plot = data['flux_B6'][ind]

F1 = np.linspace(np.min(B3_plot), np.max(B3_plot), 10)
alpha1_F2 = F1*(B6nu/B3nu)
alpha15_F2 = F1*((B6nu/B3nu)**1.5)
alpha2_F2 = F1*((B6nu/B3nu)**2)
alpha25_F2 = F1*((B6nu/B3nu)**2.5)
alpha3_F2 = F1*((B6nu/B3nu)**3)

fig1 = plt.figure()
plt.xlabel('B3 flux')
plt.ylabel('B6 flux')
plt.loglog(B3_plot, B6_plot, linestyle='', marker='*')
plt.loglog(F1, alpha1_F2, linestyle=':',label='alpha=1')
plt.loglog(F1, alpha15_F2, linestyle=':',label='alpha=1.5')
plt.loglog(F1, alpha2_F2, linestyle=':',label='alpha=2')
plt.loglog(F1, alpha25_F2, linestyle=':',label='alpha=2.5')
plt.loglog(F1, alpha3_F2, linestyle=':',label='alpha=3')
plt.legend()
#plt.show()

#B3 Kmag plot
ind = np.concatenate((np.where(np.isnan(data['flux_B3'])== False)[0],np.where(np.isnan(data['Kmag'])==False)[0]))
B3_plot_kmag = data['flux_B3'][ind]
kmag_plot_B3 = data['Kmag'][ind]

fig2 = plt.figure()
plt.xlabel('B3 flux')
plt.ylabel('K-band magnitude')
plt.semilogx(B3_plot_kmag, kmag_plot_B3, linestyle='', marker='*')
plt.ylim(15,9)

#B3 Hmag plot
ind = np.concatenate((np.where(np.isnan(data['flux_B3'])== False)[0],np.where(np.isnan(data['Hmag'])==False)[0]))
B3_plot_hmag = data['flux_B3'][ind]
hmag_plot_B3 = data['Hmag'][ind]

fig3 = plt.figure()
plt.xlabel('B3 flux')
plt.ylabel('H-band magnitude')
plt.semilogx(B3_plot_hmag, hmag_plot_B3, linestyle='', marker='*')
plt.ylim(17,9)
plt.show()
