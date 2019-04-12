from astropy.table import Table, Column
from astropy.io import fits, ascii
from astropy.modeling import blackbody
from astropy import constants
from astropy.coordinates import SkyCoord
import radio_beam
import numpy as np
import astropy.units as u


def RA_to_deg(HH, MM, SS):
            return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
    return DD + (MM/60 + SS/(60*60))*np.sign(DD)
                

FWHM_TO_SIGMA = 1/np.sqrt(8*np.log(2))

data =  Table.read('../tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')
nonconv_data = Table.read('../tables/r0.5_catalog_nonconv_apflux_final.fits')


#calculate quantities for each band
bands = ['B3', 'B6', 'B7']
imgs = ['/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
nonconv_imgs = ['/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
freqs = [98, 223.5 ,339.7672758867]


#calculated quantities not requiring loops
A = np.log10(freqs[1]) - np.log10(freqs[0])
alpha_B3B6 = (np.log10(data['ap_flux_B6'])-np.log10(data['ap_flux_B3']))/A
alpha_B3B6_err = np.sqrt((data['ap_flux_err_B6']/(A*np.log(10)*data['ap_flux_B6']))**2 + (data['ap_flux_err_B3']/(A*np.log(10)*data['ap_flux_B3']))**2)
B = np.log10(freqs[2]) - np.log10(freqs[1])
alpha_B6B7 = (np.log10(data['ap_flux_B7'])-np.log10(data['ap_flux_B6']))/B
alpha_B6B7_err = np.sqrt((data['ap_flux_err_B6']/(B*np.log(10)*data['ap_flux_B6']))**2 + (data['ap_flux_err_B7']/(B*np.log(10)*data['ap_flux_B7']))**2)

#calclated quantites for each band

#constants for mass calc
Tdust = 20*u.K
kappa0 = 2*u.cm**2/u.g
dist = 414*u.pc

tables = []
for i in range(2):
    int_flux_arrs = []
    int_flux_err_arrs = []
    lower_lum_arrs = []
    mass_arrs = []
    mass_err_arrs = []
    inclination_arrs = []
    inclination_err_arrs = []

    if i == 1:
        data = data_nonconv
    for b in range(len(bands)):
        band = bands[b]
        inclination = np.arccos(data['fwhm_min_deconv_'+band]/data['fwhm_maj_deconv_'+band])
        inclination_arrs.append((inclination*u.rad).to(u.deg).value)
        inclination_err = np.sqrt(
            (data['fwhm_min_deconv_err_'+band]**2/(data['fwhm_maj_deconv_'+band]**2 - data['fwhm_min_deconv_'+band]**2)) +
            ((data['fwhm_maj_deconv_err_'+band]**2)*(data['fwhm_min_deconv_'+band]**2)/(data['fwhm_maj_deconv_'+band]**2*(data['fwhm_maj_deconv_'+band]**2 - data['fwhm_min_deconv_'+band]**2))))
        inclination_err_arrs.append((inclination_err*u.radian).to(u.deg).value)

        fl = fits.open(imgs[b])
        fl_nonconv = fits.open(nonconv_imgs[b])
        nonconv_beam = radio_beam.Beam.from_fits_header(fl_nonconv[0].header)
        beam = radio_beam.Beam.from_fits_header(fl[0].header)

        int_flux = ((2*np.pi*data['gauss_amp_'+band]*data['fwhm_maj_'+band]*u.arcsec*data['fwhm_min_'+band]*u.arcsec*(FWHM_TO_SIGMA**2))/nonconv_beam.sr).decompose()
        int_flux_arrs.append(int_flux)
        int_flux_err = int_flux*np.sqrt((data['gauss_amp_err_'+band]/data['gauss_amp_'+band])**2+(data['fwhm_maj_err_'+band]/data['fwhm_maj_'+band])**2+(data['fwhm_min_err_'+band]/data['fwhm_min_'+band])**2)
        int_flux_err_arrs.append(int_flux_err)

        R = (np.sin(beam.minor)*dist).to(u.au)/2
        T_B = (data['ap_flux_'+band]*u.Jy).to(u.K, beam.jtok_equiv(freqs[b]*u.GHz))
        L = (4 * np.pi * R**2 * constants.sigma_sb * (T_B)**4).to(u.L_sun)
        lower_lum_arrs.append(L)

        Bnu = blackbody.blackbody_nu(freqs[b]*u.GHz, 20*u.K)
        Dmass = (data['ap_flux_'+band]*dist**2)/(kappa0*Bnu)
        Dmass_err = (data['ap_flux_err_'+band]*dist**2)/(kappa0*Bnu)
        mass_arrs.append(Dmass)
        mass_err_arrs.append(Dmass_err)

    calc_tab = Table([int_flux_arrs[0], int_flux_err_arrs[0], int_flux_arrs[1], int_flux_err_arrs[1],
                  int_flux_arrs[2], int_flux_err_arrs[2], inclination_arrs[0], inclination_err_arrs[0],
                  inclination_arrs[1], inclination_err_arrs[1], inclination_arrs[2], inclination_err_arrs[2],
                  lower_lum_arrs[0], lower_lum_arrs[1], lower_lum_arrs[2], mass_arrs[0], mass_err_arrs[0],
                  mass_arrs[1], mass_err_arrs[1], mass_arrs[2], mass_err_arrs[2], alpha_B3B6, alpha_B3B6_err,
                  alpha_B6B7, alpha_B6B7_err],
                  names=['int_flux_B3', 'int_flux_err_B3', 'int_flux_B6', 'int_flux_err_B6',
                         'int_flux_B7', 'int_flux_err_B7', 'inclination_B3', 'inclination_err_B3',
                         'inclination_B6', 'inclination_err_B6', 'inclination_B7', 'inclination_err_B7',
                         'lower_lum_B3', 'lower_lum_B6', 'lower_lum_B7', 'dust_mass_B3', 'dust_mass_err_B3',
                         'dust_mass_B6', 'dust_mass_err_B6', 'dust_mass_B7', 'dust_mass_err_B7', 'alpha_B3B6',
                         'alpha_B3B6_err', 'alpha_B6B7', 'alpha_B6B7_err'],
                 dtype=['f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8',
                        'f8','f8','f8','f8','f8','f8','f8','f8','f8'])    
    tables.append(calc_tab)



nonconv_srcs = [3,13]
for src in nonconv_srcs[0]:
    ind_tab1 = np.where(table1['D_ID'] == src)[0]
    ind_nonconv = np.where(table1['D_ID'] == src)[0]
    for col in table1.colnames:
        if col in nonconv_data.colnames:
            table1[col][ind_tab1] = nonconv_data[col][ind_nonconv]
                
    
EisnerID = Column(np.array(np.repeat('none', len(data)), dtype='S10'), name='Eisner_ID')

table1 = Table((data['D_ID'], EisnerID, data['RA_B3'], data['DEC_B3'], data['RA_err_B3'],
                data['DEC_err_B3'], data['ap_flux_B3'], data['ap_flux_err_B3'], data['ap_flux_B6'],
                data['ap_flux_err_B6'], data['ap_flux_B7'], data['ap_flux_err_B7'], data['gauss_amp_B3'],
                data['gauss_amp_err_B3'], data['gauss_amp_B6'], data['gauss_amp_err_B6'],
                data['gauss_amp_B7'], data['gauss_amp_err_B7'], calc_tab['int_flux_B3'],
                calc_tab['int_flux_err_B3'], calc_tab['int_flux_B6'], calc_tab['int_flux_err_B6'],
                calc_tab['int_flux_B7'], calc_tab['int_flux_err_B7'], data['fwhm_maj_B3'],
                data['fwhm_maj_err_B3'], data['fwhm_min_B3'],
                data['fwhm_min_err_B3'], data['pa_B3'], data['pa_err_B3'],
                data['fwhm_maj_B6'], data['fwhm_maj_err_B6'], data['fwhm_min_B6'],
                data['fwhm_min_err_B6'], data['pa_B6'], data['pa_err_B6'],
                data['fwhm_maj_B7'], data['fwhm_maj_err_B7'], data['fwhm_min_B7'],
                data['fwhm_min_err_B7'], data['pa_B7'], data['pa_err_B7'],
                data['fwhm_maj_deconv_B3'], data['fwhm_min_deconv_B3'], data['pa_deconv_B3'], 
                data['fwhm_maj_deconv_B6'], data['fwhm_min_deconv_B6'], data['pa_deconv_B6'], 
                data['fwhm_maj_deconv_B7'], data['fwhm_min_deconv_B7'], data['pa_deconv_B7'],
                calc_tab['inclination_B3'], calc_tab['inclination_err_B3'], calc_tab['inclination_B6'],
                calc_tab['inclination_err_B6'], calc_tab['inclination_B7'], calc_tab['inclination_err_B7'],
                calc_tab['alpha_B3B6'], calc_tab['alpha_B3B6_err'], calc_tab['alpha_B6B7'],
                calc_tab['alpha_B6B7_err']))

opt_depth = np.repeat('-', len(data['D_ID']))

table2 = Table((data['D_ID'], calc_tab['lower_lum_B3'], calc_tab['lower_lum_B6'], calc_tab['lower_lum_B7'],
                calc_tab['dust_mass_B3'], calc_tab['dust_mass_err_B3'], calc_tab['dust_mass_B6'],
                calc_tab['dust_mass_err_B6'], calc_tab['dust_mass_B7'], calc_tab['dust_mass_err_B7'],
                opt_depth))

B6ind = np.where(np.isnan(data['RA_B6']) == False)[0]
B7ind = np.where(np.isnan(data['RA_B7']) == False)[0]
allband_ind = np.intersect1d(B6ind, B7ind)

eisner_tab = ascii.read('../tables/eisner_tbl.txt', format='tab', delimiter='\t')
eisner_tab.remove_column('remove')

eisner_ra = []
eisner_dec = []
for row in eisner_tab:
    eisner_ra.append(RA_to_deg(float(row['alpha'][0:2]), float(row['alpha'][2:4]), float(row['alpha'][4:])))
    eisner_dec.append(DEC_to_deg(float(row['delta'][0:3]), float(row['delta'][3:5]), float(row['delta'][5:])))


    
coord_tab = SkyCoord(ra=table1['RA_B3']*u.degree, dec=table1['DEC_B3']*u.degree)

eisner_coord = SkyCoord(ra=eisner_ra*u.degree, dec=eisner_dec*u.degree)

idx, d2d, d3d = eisner_coord.match_to_catalog_sky(coord_tab)
#idx is list of indices of table1 with locations corresponding to eisner_tab
matches = np.where(d2d.value < 0.5*(1/3600))[0] #matches within 0.1 arcsec
for mat in matches:
    table1[idx[mat]]['Eisner_ID'] = eisner_tab[mat]['ID']

matched_inds = np.where(table1['Eisner_ID'] != 'none')[0]

table1 = table1[allband_ind]
table2 = table2[allband_ind]


#to do:
#change source 3 and 13 quantities with non convolved images
#fill column for whether qualitatively optically thick/thin
#rerun astrodendro to pick up on missed sources

#ascii.write(table1, '/users/jotter/summer_research_2018/tables/measured_vals_table1.fits', format='latex')
#ascii.write(table2, '/users/jotter/summer_research_2018/tables/inferred_vals_table2.fits', format='latex')
