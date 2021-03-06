from astropy.table import Table, Column
from astropy.io import fits, ascii
from astropy.modeling import blackbody
from astropy import constants
from astropy.coordinates import SkyCoord
import radio_beam
import numpy as np
import astropy.units as u


FWHM_TO_SIGMA = 1/np.sqrt(8*np.log(2))

data =  Table.read('../tables/r0.5_catalog_bgfit_may21_ulim.fits')


#calculate quantities for each band
bands = ['B3', 'B6', 'B7']
imgs = ['/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits', '/home/jotter/nrao/images/B6_convolved_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits', '/home/jotter/nrao/images/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
nonconv_imgs = ['/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits', '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits', '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
freqs = [98, 223.5 ,339.7672758867]
nonconv_srcs = [[], [30,36,45,56,63,68,70,71], [14,16,29,30,45,63,66,68,69,70]] #B3_Seq

nonconv_src_inds = []
for i in range(len(bands)):
    nonconv_ind = []
    for src in nonconv_srcs[i]:
        nonconv_ind.append(np.where(data['B3_Seq'] == src)[0][0])
    nonconv_src_inds.append(np.array(nonconv_ind))
#calclated quantites for each band

#constants for mass calc
Tdust = 20*u.K
kappa0 = 2*u.cm**2/u.g
dist = 400*u.pc
nu0 = constants.c/(1.3*u.mm)

#ind3 = np.where(nonconv_data['D_ID'] == 3)[0]
#ind13 = np.where(nonconv_data['D_ID'] == 13)[0]
#src_inds = np.concatenate((ind3, ind13))


int_flux_arrs = []
int_flux_err_arrs = []
lower_lum_arrs = []
mass_arrs = []
mass_err_arrs = []
inclination_arrs = []
inclination_err_arrs = []

sys_err = 0.20 #20% systematic error

b3flux = data['ap_flux_B3'].data
b3flux_err = np.sqrt(data['ap_flux_err_B3'].data**2 + (sys_err*b3flux)**2)
b6flux = data['ap_flux_B6'].data
b6flux_err = np.sqrt(data['ap_flux_err_B6'].data**2 + (sys_err*b6flux)**2)
b7flux = data['ap_flux_B7'].data
b7flux_err = np.sqrt(data['ap_flux_err_B7'].data**2 + (sys_err*b7flux)**2)

#calculated quantities not requiring loops
A = np.log10(freqs[1]) - np.log10(freqs[0])
alpha_B3B6 = (np.log10(b6flux)-np.log10(b3flux))/A
alpha_B3B6_err = np.sqrt((b6flux_err/(A*np.log(10)*b6flux))**2 + (b3flux_err/(A*np.log(10)*b3flux))**2)
alpha_B3B6_err2 = (1/(A*np.log(10)))*np.sqrt((b6flux_err/b6flux)**2 + (b3flux_err/b3flux)**2)

B = np.log10(freqs[2]) - np.log10(freqs[1])
alpha_B6B7 = (np.log10(b7flux)-np.log10(b6flux))/B
alpha_B6B7_err = np.sqrt((b6flux_err/(B*np.log(10)*b6flux))**2 + (b7flux_err/(B*np.log(10)*b7flux))**2)
#alpha_B3B6_err = (1/np.log(10))*np.sqrt((data['ap_flux_err_B6']/data['ap_flux_B6'])**2 + (data['ap_flux_err_B3']/data['ap_flux_B3'])**2)
    
for b in range(len(bands)):
    band = bands[b]
    inclination = np.arccos(data['fwhm_min_deconv_'+band]/data['fwhm_maj_deconv_'+band])*u.radian
    inclination_arrs.append(inclination.to(u.degree))
    inclination_err = np.sqrt((data['fwhm_min_deconv_err_'+band]**2/(data['fwhm_maj_deconv_'+band]**2- data['fwhm_min_deconv_'+band]**2)) +
                              ((data['fwhm_maj_deconv_err_'+band]**2)*(data['fwhm_min_deconv_'+band]**2) /
                               (data['fwhm_maj_deconv_'+band]**2 * (data['fwhm_maj_deconv_'+band]**2 - data['fwhm_min_deconv_'+band]**2))))*u.radian
    inclination_err_arrs.append(inclination_err.to(u.degree))

    fl = fits.open(imgs[b])
    fl_nonconv = fits.open(nonconv_imgs[b])
    nonconv_beam = radio_beam.Beam.from_fits_header(fl_nonconv[0].header)
    beam = radio_beam.Beam.from_fits_header(fl[0].header)

    int_flux = ((2*np.pi*data['gauss_amp_'+band]*data['fwhm_maj_'+band]*u.arcsec*data['fwhm_min_'+band]*u.arcsec*(FWHM_TO_SIGMA**2))/nonconv_beam.sr).decompose()
    int_flux_arrs.append(int_flux)
    int_flux_err = int_flux*np.sqrt((data['gauss_amp_err_'+band]/data['gauss_amp_'+band])**2+(data['fwhm_maj_err_'+band]/data['fwhm_maj_'+band])**2+(data['fwhm_min_err_'+band]/data['fwhm_min_'+band])**2)
    int_flux_err_arrs.append(int_flux_err)

    R = (np.sin(beam.minor)*dist).to(u.au)/2
    #T_B1 = (data['gauss_amp_'+band]*u.Jy).to(u.K, beam.jtok_equiv(freqs[b]*u.GHz))
    Fnu = data['gauss_amp_'+band]*u.Jy/beam.sr.value
    T_B = ((constants.h*(freqs[b]*u.GHz))/(constants.k_B*np.log((2*constants.h*(freqs[b]*u.GHz)**3)/(Fnu*constants.c**2).decompose() + 1))).decompose()

    L = (4 * np.pi * R**2 * constants.sigma_sb * (T_B)**4).to(u.L_sun)

    R_alt = (np.sin(nonconv_beam.minor)*dist).to(u.au)/2
    Fnu_alt = data['gauss_amp_'+band][nonconv_src_inds[i]]*u.Jy/nonconv_beam.sr.value
    T_B_alt = ((constants.h*(freqs[b]*u.GHz))/(constants.k_B*np.log((2*constants.h*(freqs[b]*u.GHz)**3)/(Fnu_alt*constants.c**2).decompose() + 1))).decompose()
    L_alt = (4 * np.pi * R_alt**2 * constants.sigma_sb * (T_B_alt)**4).to(u.L_sun)

    #print(L[nonconv_src_inds[i]])
    
    L[nonconv_src_inds[i]] = L_alt

    #print(L[nonconv_src_inds[i]])
    #print(L)
    lower_lum_arrs.append(L)

    
    
    Bnu = blackbody.blackbody_nu(freqs[b]*u.GHz, 20*u.K)
    Dmass = (data['ap_flux_'+band]*u.Jy*dist**2)/(kappa0*(freqs[b]*u.GHz/nu0)*Bnu)
    Dmass_err = (data['ap_flux_err_'+band]*u.Jy*dist**2)/(kappa0*(freqs[b]*u.GHz/nu0)*Bnu)
    Dmass = (Dmass.decompose()).to(u.earthMass*u.sr)
    Dmass_err = (Dmass_err.decompose()).to(u.earthMass*u.sr)
        
    mass_arrs.append(Dmass.value)
    mass_err_arrs.append(Dmass_err.value)

Bnu_B6 = blackbody.blackbody_nu(freqs[1]*u.GHz, 20*u.K)
B6_dmass_ulim = ((data['B6_flux_ulim'].data*u.Jy*dist**2)/(kappa0*(freqs[1]*u.GHz/nu0)*Bnu)).to(u.earthMass*u.sr).value
Bnu_B7 = blackbody.blackbody_nu(freqs[2]*u.GHz, 20*u.K)
B7_dmass_ulim = ((data['B7_flux_ulim'].data*u.Jy*dist**2)/(kappa0*(freqs[2]*u.GHz/nu0)*Bnu)).to(u.earthMass*u.sr).value

tab = Table([data['ID'].data, int_flux_arrs[0], int_flux_err_arrs[0], int_flux_arrs[1],
             int_flux_err_arrs[1], int_flux_arrs[2], int_flux_err_arrs[2], inclination_arrs[0],
             inclination_err_arrs[0], inclination_arrs[1], inclination_err_arrs[1], inclination_arrs[2],
             inclination_err_arrs[2], lower_lum_arrs[0], lower_lum_arrs[1], lower_lum_arrs[2],
             mass_arrs[0], mass_err_arrs[0], mass_arrs[1], mass_err_arrs[1], mass_arrs[2],
             mass_err_arrs[2], alpha_B3B6, alpha_B3B6_err, alpha_B6B7, alpha_B6B7_err, B6_dmass_ulim, B7_dmass_ulim],
             names=['ID', 'int_flux_B3', 'int_flux_err_B3', 'int_flux_B6', 'int_flux_err_B6',
                    'int_flux_B7', 'int_flux_err_B7', 'inclination_B3', 'inclination_err_B3',
                    'inclination_B6', 'inclination_err_B6', 'inclination_B7', 'inclination_err_B7',
                    'lower_lum_B3', 'lower_lum_B6', 'lower_lum_B7', 'dust_mass_B3', 'dust_mass_err_B3',
                    'dust_mass_B6', 'dust_mass_err_B6', 'dust_mass_B7', 'dust_mass_err_B7', 'alpha_B3B6',
                    'alpha_B3B6_err', 'alpha_B6B7', 'alpha_B6B7_err','dust_mass_ulim_B6', 'dust_mass_ulim_B7'],
             dtype=['i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8',
                    'f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'])    

tab.write('/home/jotter/nrao/summer_research_2018/tables/r0.5_may21_calc_vals.fits', format='fits', overwrite=True)
