from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt
import numpy as np

def disk_size_hist(arrs, labels, filename):
    fig = plt.figure()
    for a in range(len(arrs)):
        size_arr = arrs[a]
        size_arr = size_arr[np.isnan(size_arr)==False]
        print(len(size_arr))
        hist, bins = np.histogram(size_arr, density=False)
        plotpts = []
        widths = []
        for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
            plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
            widths.append((bins[b+1]-bins[b])/2)
    
        plt.bar(plotpts, hist, widths, edgecolor = 'black', label=labels[a])
        plt.legend()
        plt.xlabel('fwhm major (as)')
        plt.ylabel('number of disks')
        plt.savefig('plots/size_plots'+filename, dpi=500)
        
def img_size_comp(arrs, err_arrs,  labels):
    fig = plt.figure()
    inds1 = np.where(err_arrs[0] < 0.2)[0]
    inds2 = np.where(err_arrs[1] < 0.2)[0]
    inds = np.intersect1d(inds1, inds2)
    plt.errorbar(arrs[0][inds], arrs[1][inds], xerr=err_arrs[0][inds], yerr=err_arrs[1][inds], linestyle='', marker='o')
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
