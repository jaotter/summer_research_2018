from astropy.table import Table
from astropy.io import ascii
from scipy.stats import kstest
import astropy.units as u

import math
import matplotlib.pyplot as plt
import numpy as np

dataB3 = Table.read('../tables/table_meas_B3.fits')

ind = np.where(np.isnan(dataB3['inclination_B3']) == False)[0]
incl_bins = np.arange(0,91,10)
hist, bins = np.histogram(dataB3['inclination_B3'][ind], bins=incl_bins)

plotpts = []
widths = []
prob = []
for b in range(len(incl_bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(incl_bins[b] + (incl_bins[b+1]-incl_bins[b])/2)
    widths.append((incl_bins[b+1]-incl_bins[b]))

#add sin(i) histogram
tot = np.sum(hist)
prob = np.sin(plotpts*u.degree)
prob = prob/np.sum(prob)
sin_hist = tot*prob


incl = dataB3['inclination_B3'][ind]*u.degree
incl = incl.to(u.radian)

Dval, pval = kstest(incl, np.sin)
print('KS statistic: %f, p value: %f' % (Dval, pval))

Dval, pval = kstest(incl, 'cosine')
print('KS statistic: %f, p value: %f' % (Dval, pval))


plt.figure()
plt.bar(plotpts, sin_hist, widths, edgecolor='k', alpha=.5, label=r'$\sin(i)$')
plt.bar(plotpts, hist, widths, edgecolor='k', alpha=0.5, label='band 3')
plt.xlabel('Band 3 inclination angle (degrees)')
plt.ylabel('number')
plt.legend()
plt.xlim(0,90)
plt.savefig('plots/incl_hist.png',dpi=300)


#histogram of cos(i)

cos_bins = np.linspace(0,1,10)
cos_incl_hist, bins = np.histogram(np.cos(dataB3['inclination_B3'][ind]*u.degree), bins=cos_bins)

plotpts = []
widths = []
for b in range(len(cos_bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(cos_bins[b] + (cos_bins[b+1]-cos_bins[b])/2)
    widths.append((cos_bins[b+1]-cos_bins[b]))

#add uniform histogram
tot = np.sum(cos_incl_hist)
unif_hist= np.repeat(tot/len(cos_incl_hist), len(cos_incl_hist))

Dval, pval = kstest(np.cos(incl), 'uniform')
print('KS statistic: %f, p value: %f' % (Dval, pval))


plt.figure()
plt.bar(plotpts, unif_hist, widths, edgecolor='k', alpha=.5, label='uniform')
plt.bar(plotpts, cos_incl_hist, widths, edgecolor='k', alpha=0.5, label='band 3')
plt.xlabel('cos(i)')
plt.ylabel('number')
plt.legend()
plt.xlim(0,1)
plt.savefig('plots/cos_incl_hist.png',dpi=300)



