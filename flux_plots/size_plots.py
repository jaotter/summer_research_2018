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
            widths.append((bins[b+1]-bins[b]))
    
        plt.bar(plotpts, hist, widths, edgecolor = 'black', label=labels[a], alpha=0.5)
    plt.legend()
    plt.xlabel('deconvolved fwhm major (as)')
    plt.ylabel('number of disks')
    plt.xlim(0,0.35)
    plt.savefig('plots/size_plots/'+filename, dpi=500)
        
def size_comp(conv_arrs, deconv_arrs, conv_errs, deconv_errs, labels, filename):
    #plot sizes of disks in two bands, arr1 and arr2 should have same length and sources
    #labels - arr1 label, then arr2

    fig = plt.figure()

    deconv_inds = np.where(np.isnan(deconv_arrs[0]) == False)[0]
    for ind in deconv_inds:
        plt.errorbar([deconv_arrs[0][ind], conv_arrs[0][ind]], [deconv_arrs[1][ind], conv_arrs[1][ind]], xerr=[deconv_errs[0][ind], conv_errs[0][ind]], yerr=[deconv_errs[1][ind], conv_errs[1][ind]], marker='o', color='r')

    plt.errorbar(conv_arrs[0], conv_arrs[1], xerr=conv_errs[0], yerr=conv_errs[1], marker='o',linestyle='', label='convolved sizes')
    
        
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.xlim(0.0,0.3)
    plt.ylim(0.0,0.3)
    plt.plot(np.arange(0,1,0.1), np.arange(0,1,0.1), color='k')
    plt.legend()
    plt.savefig('plots/size_plots/'+filename, dpi=400)

def size_comp_simple(arrs, errs, labels, filename):
    
    #plot sizes of disks in two bands, arr1 and arr2 should have same length and sources
    #labels - arr1 label, then arr2

    fig = plt.figure()

    plt.errorbar(arrs[0], arrs[1], xerr=errs[0], yerr=errs[1], linestyle='', marker='.')    
        
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.xlim(0.0,0.3)
    plt.ylim(0.0,0.3)
    plt.plot(np.arange(0,1,0.1), np.arange(0,1,0.1), color='k')
    plt.savefig('plots/size_plots/'+filename, dpi=400)


    
def img_size_comp(arrs, err_arrs,  labels):
    fig = plt.figure()
    inds1 = np.where(err_arrs[0] < 0.2)[0]
    inds2 = np.where(err_arrs[1] < 0.2)[0]
    inds = np.intersect1d(inds1, inds2)
    plt.errorbar(arrs[0][inds], arrs[1][inds], xerr=err_arrs[0][inds], yerr=err_arrs[1][inds], linestyle='', marker='o')
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
