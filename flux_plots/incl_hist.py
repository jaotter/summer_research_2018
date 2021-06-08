from astropy.table import Table
from astropy.io import ascii
from scipy.stats import ks_1samp, ks_2samp, gaussian_kde, uniform
import astropy.units as u

import math
import matplotlib.pyplot as plt
import numpy as np

def f(x): #for 2 sample ks test
    return (1 - np.cos(x))*(x>0)*(x<np.pi/2)+(x>np.pi/2)


dataB3 = Table.read('../tables/r0.5_may21_calc_vals_mask.fits')
#dataB3 = Table.read('../tables/r0.5_jun20_calc_vals.fits')

ind = np.where(np.isnan(dataB3['inclination_B3']) == False)[0]

data = Table.read('../tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
#incl = data['fwhm_maj_deconv_B3'] / data['fwhm_min_deconv_B3']

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


incl = dataB3['inclination_B3'][ind]
incl_rad = incl.to(u.radian).value
incl_filtered = incl_rad % np.pi/2

x = np.random.rand(10000)
y = np.arccos(1-x)

#2 sample ks test
Dval, pval = ks_2samp(incl_rad, y) 
print('2 sample KS statistic: %f, p value: %f' % (Dval, pval))

#1 sample test
Dval, pval = ks_1samp(np.cos(incl_rad), uniform.cdf) 
print('1 sample KS statistic: %f, p value: %f' % (Dval, pval))

#1 sample test
Dval, pval = ks_1samp(incl_rad, lambda x: 1-np.cos(x)) 
print('1 sample KS statistic: %f, p value: %f' % (Dval, pval))


minor = data['fwhm_min_deconv_B3'][np.isnan(data['fwhm_min_deconv_B3'])==False]
major = data['fwhm_maj_deconv_B3'][np.isnan(data['fwhm_min_deconv_B3'])==False]
inclination = np.arccos(minor/major)
print(ks_1samp(np.cos(inclination), uniform.cdf))


#KDE
#cos_grid = np.linspace(0,1,100)
degree_grid = np.linspace(0,90,100)
#cos_kde = gaussian_kde(np.cos(dataB3['inclination_B3'][ind]*u.degree))
degree_kde = gaussian_kde(dataB3['inclination_B3'][ind], bw_method='scott')
#cos_pdf = cos_kde.evaluate(cos_grid)
#norm_cos_pdf = cos_pdf*len(ind)*(1/9)
#to normalize, multiply by bin width and number of points
degree_pdf = degree_kde.evaluate(degree_grid)
norm_degree_pdf = degree_pdf*widths[0]*np.sum(hist)

plt.figure()
plt.bar(plotpts, sin_hist, widths, edgecolor='k', alpha=.5, label=r'$\sin(i)$')
plt.bar(plotpts, hist, widths, edgecolor='k', alpha=0.5, label='band 3')
plt.plot(degree_grid, norm_degree_pdf, color='tab:red', label='KDE')
plt.xlabel('Band 3 inclination angle (degrees)')
plt.ylabel('number')
plt.legend()
plt.xlim(0,90)
plt.savefig('../plots/incl_hist.pdf',dpi=300)

# attempt to estimate scale height
# as a function of scale height, what is the p-value?  What p-value would imply a 2-sigma consistency?  3-sigma?
scaleheights = np.linspace(0,10,100)
b3sh = []
distance_pc = 400
for scaleheight in scaleheights:
    inclination = np.arccos((minor**2 - (scaleheight/distance_pc)**2)**0.5/major)
    inclination = inclination[np.isfinite(inclination)]
    D,P = (ks_1samp(inclination, lambda x: 1-np.cos(x)))
    b3sh.append(P)
plt.figure()
plt.plot(scaleheights, b3sh)
plt.hlines([0.05,0.003], 0, 10)
plt.xlabel("Scaleheight")
plt.ylabel("P-value")
# the result is kind of marginal - between 2 and 3 sigma




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
