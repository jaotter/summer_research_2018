from astropy.table import Table
from astropy.io import ascii
from scipy.stats import kstest, ks_2samp, gaussian_kde
import astropy.units as u

import math
import matplotlib.pyplot as plt
import numpy as np

def f(x): #for 2 sample ks test
    return (1 - np.cos(x))*(x>0)*(x<np.pi/2)+(x>np.pi/2)


dataB3 = Table.read('../tables/r0.5_mar21_calc_vals.fits')
#dataB3 = Table.read('../tables/r0.5_jun20_calc_vals.fits')

#ind = np.where(np.isnan(dataB3['inclination_B3']) == False)[0]
ind = np.where(dataB3['inclination_B3'] < 1000)[0]


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
incl_rad = incl.to(u.radian).value
incl_filtered = incl_rad % np.pi/2

x = np.random.rand(10000)
y = np.arccos(1-x)

#2 sample ks test
Dval, pval = ks_2samp(incl_filtered, y) 
print('2 sample KS statistic: %f, p value: %f' % (Dval, pval))

#1 sample test
Dval, pval = kstest(incl_filtered, f) 
print('1 sample KS statistic: %f, p value: %f' % (Dval, pval))

print(dataB3['inclination_B3'][ind])

#KDE
#cos_grid = np.linspace(0,1,100)
degree_grid = np.linspace(0,90,100)
#cos_kde = gaussian_kde(np.cos(dataB3['inclination_B3'][ind]*u.degree))
degree_kde = gaussian_kde(dataB3['inclination_B3'][ind], bw_method='scott')
#cos_pdf = cos_kde.evaluate(cos_grid)
#norm_cos_pdf = cos_pdf*len(ind)*(1/9)
#to normalize, multiply by bin width and number of points
degree_pdf = degree_kde.evaluate(degree_grid)
print(widths[0])
norm_degree_pdf = degree_pdf*widths[0]*np.sum(hist)
print(norm_degree_pdf)

plt.figure()
plt.bar(plotpts, sin_hist, widths, edgecolor='k', alpha=.5, label=r'$\sin(i)$')
plt.bar(plotpts, hist, widths, edgecolor='k', alpha=0.5, label='band 3')
plt.plot(degree_grid, norm_degree_pdf, color='tab:red', label='KDE')
plt.xlabel('Band 3 inclination angle (degrees)')
plt.ylabel('number')
plt.legend()
plt.xlim(0,90)
plt.savefig('/home/jotter/nrao/plots/incl_hist.pdf',dpi=300)



#histogram of cos(i)
'''
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

Dval, pval = kstest(np.cos(incl_filtered), 'uniform')
print('KS statistic: %f, p value: %f' % (Dval, pval))



plt.figure()
plt.bar(plotpts, unif_hist, widths, edgecolor='k', alpha=.5, label='uniform')
plt.bar(plotpts, cos_incl_hist, widths, edgecolor='k', alpha=0.5, label='band 3')
plt.plot(cos_grid, norm_cos_pdf, color='b', label='KDE')
plt.xlabel('cos(i)')
plt.ylabel('number')
plt.legend()
plt.xlim(0,1)
plt.savefig('../../plots/cos_incl_hist.png',dpi=300)



'''
