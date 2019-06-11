from astropy.table import Table
from astropy.io import ascii

import matplotlib.pyplot as plt
import numpy as np

dataB3 = Table.read('../tables/table_meas_B3.fits')

ind = np.where(np.isnan(dataB3['inclination_B3']) == False)[0]
hist, bins = np.histogram(dataB3['inclination_B3'][ind], bins=9)

plotpts = []
widths = []
for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))

#add cos(i) histogram
tot = np.sum(hist)
prob = np.abs(np.cos(plotpts))
prob = prob/np.sum(prob)
cos_hist = tot*prob

print(np.sum(cos_hist), np.sum(hist))

plt.bar(plotpts, cos_hist, widths, edgecolor='k', alpha=.5, label=r'$\cos(i)$')


plt.bar(plotpts, hist, widths, edgecolor='k', alpha=0.5, label='band 3')
plt.xlabel('Band 3 inclination angle (degrees)')
plt.ylabel('number')
plt.legend()
plt.savefig('plots/incl_hist.png',dpi=300)
