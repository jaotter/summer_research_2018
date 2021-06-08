from astropy.io import fits
from radio_beam import Beam
from astropy.convolution import convolve_fft
from astropy.wcs import WCS
from astropy.table import Table
from matplotlib import patches
from astropy.coordinates import SkyCoord, Angle
from scipy.optimize import curve_fit
from gaussfit_catalog import gaussfit_catalog
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import regions
import sys
import os

def draw_circle(r,size,flux): #draw circle radius r in 2d rectangular array arr
    xx, yy = np.mgrid[:size, :size]
    radius = np.sqrt((xx - size/2)**2 + (yy - size/2)**2)
    circle_arr = np.zeros((size, size))
    circle_ind = np.where(radius < r)
    circle_arr[circle_ind] = flux/len(circle_ind)
    return circle_arr

def f(x, a, b):
    return a*(x**b)


def size_ulim(band, radii_au, n_rep=10):
    radii_au = sorted(radii_au)[::-1]
    radii_as = ((radii_au*u.AU).to(u.pc)/(400*u.pc)*u.rad).to(u.arcsecond).value
    
    B3file = '/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits'
    B6file = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'#B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits'

    B7file = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'#B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'

    if band == 'B3':
        img = B3file
    if band == 'B6':
        img = B6file
    if band == 'B7':
        img = B7file

    orig_fl = fits.open(img)

    wcs = WCS(orig_fl[0].header).celestial
    beam = Beam.from_fits_header(orig_fl[0].header)
    
    pixel_scale = np.abs(wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    kernel = beam.as_kernel(pixel_scale)

    size_as=3*u.arcsecond #square size of fake image
    size = (size_as/pixel_scale).decompose().value
    size = int(size)
    print(size)
    
    noise_nonconv = 1e-4
    #snrs = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22.5,25,27.5,30]
    snrs = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22.5,25,27.5,30]
    #snrs = [5]
    full_upper_lim = []
    full_snrs = []
    
    for n in range(n_rep):

        upper_size_lim = []
        snr_vals = []
        for id in range(len(snrs)):
            snr = snrs[id] + np.random.normal(0,0.2,1)
            snr_vals.append(snr)

            print(f'Fitting S/N {snr} disk')

            src_name = 'SNR_'+str(snrs[id])+'_'+band

            recovered_fwhm_maj = []
           
            for i,rad_init in enumerate(radii_as):
                noise_arr_raw = np.random.normal(0, noise_nonconv, (size,size))
                noisy_arr_conv = convolve_fft(noise_arr_raw, kernel)

                noise_conv = np.std(noisy_arr_conv)

                flux = noise_conv*snr

                rand_val_au = np.random.normal(0,0.2,1)
                rand_val_as = (((rand_val_au*u.AU)/(400*u.pc)).decompose()*u.radian).to(u.arcsecond).value
                rad = rad_init + rand_val_as
                rad_au = ((rad*u.arcsecond).to(u.radian) * (400*u.pc).to(u.AU)).value
                
                print(f'Disk radius {rad} as')
                #create fake image data - uniform disk with radius r
                r_arcsec = rad*u.arcsec #radius in arcsec
                r = (r_arcsec/pixel_scale).decompose() #radius in pixel
                fake_data = draw_circle(r,size,1)
                fake_data_conv = convolve_fft(fake_data, kernel)

                disk_conv_peak = np.max(fake_data_conv)
                fake_data_conv_scaled = fake_data_conv * (flux/disk_conv_peak)
                
                noisy_conv_fake_img = fake_data_conv_scaled + noisy_arr_conv

                print(f'std raw noise: {np.nanstd(noise_arr_raw)}, std conv noise: {np.nanstd(noisy_arr_conv)}, std image: {np.nanstd(noisy_conv_fake_img)}')
                print(f'disk raw flux total: {np.sum(fake_data)}, disk conv flux total: {np.sum(fake_data_conv)}, scaled disk flux total: {np.sum(fake_data_conv_scaled)}')
                print(f'disk raw flux peak: {np.max(fake_data)}, disk conv flux peak: {np.max(fake_data_conv)}, scaled disk flux peak: {np.max(fake_data_conv_scaled)}, theoretical flux: {flux}')


                
                savepth = f'/home/jotter/nrao/plots/size_lims/imgs/'
                if not os.path.exists(savepth):
                    os.mkdir(savepth)

                new_fl = orig_fl
                new_fl[0].data = noisy_conv_fake_img
                new_fl.writeto(f'/home/jotter/nrao/plots/size_lims/imgs/fake_conv_img_r{str(radii_au[i])}_{src_name}.fits', overwrite=True)
                
                loc = wcs.wcs_pix2world(size/2,size/2,1)
                reg = regions.CircleSkyRegion(center=SkyCoord(loc[0], loc[1], unit='deg'), radius=0.5*u.arcsec, meta={'text':f'r={str(radii_au[i])}_{n}'})
                reg_pix = reg.to_pixel(wcs)


                gaussfit = gaussfit_catalog(f'/home/jotter/nrao/plots/size_lims/imgs/fake_conv_img_r{str(radii_au[i])}_{src_name}.fits', [reg], Angle(0.3, 'arcsecond'), savepath=savepth, max_radius_in_beams=15)

                source_size = Beam(major=gaussfit[f'r={str(radii_au[i])}_{n}']['fwhm_major'], minor=gaussfit[f'r={str(radii_au[i])}_{n}']['fwhm_minor'], pa=(gaussfit[f'r={str(radii_au[i])}_{n}']['pa'].value+90)*u.degree)
                try:
                    deconv_size = source_size.deconvolve(beam)
                    recovered_fwhm_maj.append(deconv_size.major.value)
                except ValueError:
                    print('could not be deconvolved')
                    recovered_fwhm_maj.append(0)
                    upper_size_lim.append(rad_au[0])
                    break

                if i == len(radii_as)-1:
                    upper_size_lim.append(0)

                
        full_upper_lim.append(np.array(upper_size_lim).flatten())
        full_snrs.append(np.array(snr_vals).flatten())

    full_upper_lim = np.concatenate(full_upper_lim).flatten()
    full_snrs = np.concatenate(full_snrs).flatten()    
    
    params, params_cov = curve_fit(f, full_snrs, full_upper_lim)

    plt.figure()
    plt.plot(full_snrs, full_upper_lim, marker='o', linestyle='', markersize=4)
    plt.plot(snrs, f(snrs, params[0], params[1]))
    plt.xlabel('SNR')
    #plt.ylim(0,20)
    plt.ylabel('Greatest deconvolvable radius (AU)')
    plt.savefig(f'/home/jotter/nrao/plots/size_lims/SNR_radius_plot_{band}_{n_rep}.png', bbox_inches='tight')
    plt.close()

    plt.figure()
    plt.hist2d(full_snrs, full_upper_lim, cmap='Blues')
    plt.plot(snrs, f(snrs, params[0], params[1]), color='tab:orange')
    plt.xlabel('SNR')
    #plt.ylim(0,20)
    plt.ylabel('Greatest deconvolvable radius (AU)')
    plt.colorbar()
    plt.savefig(f'/home/jotter/nrao/plots/size_lims/SNR_radius_hist2d_{band}_{n_rep}.png', bbox_inches='tight')
    plt.close()

    tab = Table((full_snrs, full_upper_lim), names=('SNR','R_UL'))
    tab.write(f'/home/jotter/nrao/plots/size_lims/size_lim_table_n={n_rep}_{band}.fits', overwrite=True)

    return params, params_cov


def table_fit(table_path, band, n_rep=10, plot=True):
    tab = Table.read(table_path)
    snrs = tab['SNR']
    upper_lims = tab['R_UL']

    params, params_cov = curve_fit(f, snrs, upper_lims)

    print('fit parameters:', params, params_cov)

    if plot == True:
        plot_snrs = np.linspace(np.min(snrs), np.max(snrs), 30)
        plt.figure()
        plt.hist2d(snrs, upper_lims, cmap='Blues')
        plt.plot(plot_snrs, f(plot_snrs, params[0], params[1]), color='tab:orange', label=f'$R = {np.round(params[0],1)}\sigma^{{{np.round(params[1],1)}}}$')
        plt.xlabel('S/N')
        #plt.ylim(0,20)
        plt.ylabel('Smallest deconvolvable FWHM (AU)')
        plt.colorbar()
        plt.legend()
        plt.savefig(f'/home/jotter/nrao/plots/size_lims/SNR_radius_hist2d_{band}_{n_rep}_2.png', bbox_inches='tight')
        plt.close()

    
    return params,params_cov


def add_upper_lims(params, band, tab=None):
    if tab is None:
        tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21.fits')
    band_ind = np.where(np.isnan(tab[f'SNR_{band}']) == False)[0]

    nondeconv_ind = np.where(np.isnan(tab[f'fwhm_maj_deconv_{band}']) == True)[0]
    nd_band_ind = np.intersect1d(band_ind, nondeconv_ind)
    
    upper_lim_arr = np.repeat(np.nan, len(tab))

    print(band)
    
    for ind in nd_band_ind:
        snr = tab[f'SNR_{band}'][ind]
        ulim = f(snr, params[0], params[1])
        upper_lim_arr[ind] = ulim

        
    if f'upper_lim_{band}' in tab.columns:
        tab[f'upper_lim_{band}'] = upper_lim_arr*u.AU
    else:
        tab.add_column(upper_lim_arr*u.AU, name=f'upper_lim_{band}')

    return tab
        
#radii_au = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,35]
#radii_au_B6 = [2,3.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30]
#radii_au_B7 = [2,3.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30]

#radii_au = [1,5,10,15,20]
#params_B3, params_cov = size_ulim('B3', radii_au, n_rep=20)
#params_B6, params_cov = size_ulim('B6', radii_au_B6, n_rep=20)
#params_B7, params_cov = size_ulim('B7', radii_au_B7, n_rep=20)

#print('B3 params', params_B3)
#print('B6 params', params_B6)
#print('B7 params', params_B7)

#print('B3')
#params_B3, params_cov = table_fit('/home/jotter/nrao/plots/size_lims/size_lim_table_n=20_B3.fits', 'B3', plot=True, n_rep=20)
#print('B6')
#params_B6, params_cov = table_fit('/home/jotter/nrao/plots/size_lims/size_lim_table_n=20_B6.fits', 'B6', plot=True, n_rep=20)
#print('B7')
#params_B7, params_cov = table_fit('/home/jotter/nrao/plots/size_lims/size_lim_table_n=20_B7.fits', 'B7', plot=True, n_rep=20)


params_B3 = [48, -0.5]
tab_b3ulim = add_upper_lims(params_B3, 'B3', tab=Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits'))
params_B6 = [26.09342253, -0.57997014]
tab_b3b6ulim = add_upper_lims(params_B6, 'B6', tab=tab_b3ulim)
params_B7 = [18.32718856, -0.73984847] 
tab_ulim = add_upper_lims(params_B7, 'B7', tab=tab_b3b6ulim)

tab_ulim.write('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits', overwrite=True)

