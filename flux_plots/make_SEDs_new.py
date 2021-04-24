from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import radio_beam
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from deconv_sources import deconv_srcs
import regions
from functools import reduce
import fnmatch
from astropy.coordinates import Angle, SkyCoord

def get_ind(names): #returns indices where sources are detected in all bands
    data = fits.getdata('tables/r0.5_catalog_bgfit_apr20.fits')
    flux_inds = [np.where(np.isnan(data['ap_flux_'+n]) == False)[0] for n in names]
    #fit_ind =  np.where(data['good_fit_flag'] == True)[0]
    ind1 = reduce(np.intersect1d, flux_inds)
    #ind = np.intersect1d(ind1, fit_ind)
    ind = ind1
    return ind


def plot_SEDs(names):
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7':339.7672758867*u.GHz}
    #freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7_hr':339.7672758867*u.GHz}
    data = fits.getdata('tables/r0.5_catalog_bgfit_apr20.fits')
    ind = get_ind(names)

    alpha = [1.5,2,2.5]

    imgs = ['images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', 'images/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', 'images/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
    freq_x = np.array([freqs[n].value for n in names])
    for i in ind:
        #fluxes = [data['gauss_amp_'+n][i] for n in names]
        #flux_err = [data['gauss_amp_err_'+n][i] for n in names]

        fluxes = [data['ap_flux_'+n][i] for n in names]
        flux_err = [data['ap_flux_err_'+n][i] for n in names]

        F2 = fluxes[1]
        nu2 = freq_x[1]
        if data['D_ID'][i] == 13:
            print(fluxes)
        
        plt.figure()
        plt.errorbar(freq_x, fluxes, yerr=flux_err, linestyle='', marker='o')
               
        for j in range(len(alpha)):
            F1 = F2*(freq_x/nu2)**alpha[j]
            plt.plot(freq_x, F1, linestyle='-', label=r'$\alpha = $'+str(alpha[j]))

        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('Aperture flux (Jy)')
        #plt.ylabel('gaussian amplitude (Jy)')
        plt.xlabel('Frequency (GHz)')
        plt.legend()

        locs_x = [0.3, 0.5, 0.7]
        locs_y = [0.15, 0.15, 0.15]
        
        for j in range(len(imgs)):
            img_file = imgs[j]
            img_fl = fits.open(img_file)
            img_data = img_fl[0].data.squeeze()
            img_header = img_fl[0].header
            mywcs = WCS(img_header).celestial
            pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
            center_coord = SkyCoord(data['RA_'+names[j]][i], data['DEC_'+names[j]][i], frame='icrs', unit=(u.deg, u.deg))
            center_coord_pix = center_coord.to_pixel(mywcs)
            if j == 0:
                pix_major_fwhm = ((data['fwhm_maj_B3'][j]*u.arcsec).to(u.degree)/pixel_scale).decompose()
                size = 3*pix_major_fwhm.value
            cutout = Cutout2D(img_data, center_coord_pix, size, mywcs, mode='partial')
            a = plt.axes([locs_x[j], locs_y[j], .2, .2])
            plt.imshow(cutout.data, origin='lower')
            plt.title(names[j])
            a.tick_params(labelleft=False, labelbottom=False)

        plt.savefig('plots/SEDs/SED_'+str(data['D_ID'][i])+'_apflux_bgfit_apr20.png', dpi=300)
        plt.close()
plot_SEDs(['B3', 'B6', 'B7'])
