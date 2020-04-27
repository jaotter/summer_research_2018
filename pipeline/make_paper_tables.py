from astropy.table import Table, Column
from astropy.io import fits, ascii
from astropy.modeling import blackbody
from astropy import constants
from astropy.coordinates import SkyCoord
import radio_beam
import numpy as np
import astropy.units as u

from latex_info import *

def RA_to_deg(HH, MM, SS):
            return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
    return DD + (MM/60 + SS/(60*60))*np.sign(DD)
                

FWHM_TO_SIGMA = 1/np.sqrt(8*np.log(2))

conv_data =  Table.read('../tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')
print(len(conv_data))
print(conv_data['D_ID'][-1])
nonconv_data = Table.read('../tables/r0.5_catalog_nonconv_apflux_final.fits')


#calculate quantities for each band
bands = ['B3', 'B6', 'B7']
imgs = ['/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
nonconv_imgs = ['/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/OrionB6/Orion_SourceI_B6_continuum_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/OrionB7/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']
freqs = [98, 223.5 ,339.7672758867]


#calclated quantites for each band

#constants for mass calc
Tdust = 20*u.K
kappa0 = 2*u.cm**2/u.g
dist = 414*u.pc
nu0 = constants.c/(1.3*u.mm)

#ind3 = np.where(nonconv_data['D_ID'] == 3)[0]
#ind13 = np.where(nonconv_data['D_ID'] == 13)[0]
#src_inds = np.concatenate((ind3, ind13))

calc_tables = []

data = conv_data
for i in range(2):
    int_flux_arrs = []
    int_flux_err_arrs = []
    lower_lum_arrs = []
    mass_arrs = []
    mass_err_arrs = []
    inclination_arrs = []
    inclination_err_arrs = []

    #if i == 1:
    #    data = nonconv_data[src_inds]

    #calculated quantities not requiring loops
    A = np.log10(freqs[1]) - np.log10(freqs[0])
    alpha_B3B6 = (np.log10(data['ap_flux_B6'])-np.log10(data['ap_flux_B3']))/A
    alpha_B3B6_err = np.sqrt((data['ap_flux_err_B6']/(A*np.log(10)*data['ap_flux_B6']))**2 + (data['ap_flux_err_B3']/(A*np.log(10)*data['ap_flux_B3']))**2)
    B = np.log10(freqs[2]) - np.log10(freqs[1])
    alpha_B6B7 = (np.log10(data['ap_flux_B7'])-np.log10(data['ap_flux_B6']))/B
    alpha_B6B7_err = np.sqrt((data['ap_flux_err_B6']/(B*np.log(10)*data['ap_flux_B6']))**2 + (data['ap_flux_err_B7']/(B*np.log(10)*data['ap_flux_B7']))**2)

    
    for b in range(len(bands)):
        band = bands[b]
        inclination = np.arccos(data['fwhm_min_deconv_'+band]/data['fwhm_maj_deconv_'+band])
        inclination_arrs.append((inclination*u.rad).to(u.deg).value)
        inclination_err = np.sqrt((data['fwhm_min_deconv_err_'+band]**2/(data['fwhm_maj_deconv_'+band]**2- data['fwhm_min_deconv_'+band]**2)) +
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
        T_B1 = (data['gauss_amp_'+band]*u.Jy).to(u.K, beam.jtok_equiv(freqs[b]*u.GHz))
        Fnu = data['gauss_amp_'+band]*u.Jy/beam.sr.value
        T_B = ((constants.h*(freqs[b]*u.GHz))/(constants.k_B*np.log((2*constants.h*(freqs[b]*u.GHz)**3)/(Fnu*constants.c**2).decompose() + 1))).decompose()

        print(T_B1, T_B)
        L = (4 * np.pi * R**2 * constants.sigma_sb * (T_B)**4).to(u.L_sun)
        lower_lum_arrs.append(L)
        
        Bnu = blackbody.blackbody_nu(freqs[b]*u.GHz, 20*u.K)
        Dmass = (data['ap_flux_'+band]*u.Jy*dist**2)/(kappa0*(freqs[b]*u.GHz/nu0)*Bnu)
        Dmass_err = (data['ap_flux_err_'+band]*u.Jy*dist**2)/(kappa0*(freqs[b]*u.GHz/nu0)*Bnu)
        Dmass = (Dmass.decompose()).to(u.earthMass*u.sr)
        Dmass_err = (Dmass_err.decompose()).to(u.earthMass*u.sr)
        
        mass_arrs.append(Dmass.value)
        mass_err_arrs.append(Dmass_err.value)

    tab = Table([data['D_ID'].data, int_flux_arrs[0], int_flux_err_arrs[0], int_flux_arrs[1],
                 int_flux_err_arrs[1], int_flux_arrs[2], int_flux_err_arrs[2], inclination_arrs[0],
                 inclination_err_arrs[0], inclination_arrs[1], inclination_err_arrs[1], inclination_arrs[2],
                 inclination_err_arrs[2], lower_lum_arrs[0], lower_lum_arrs[1], lower_lum_arrs[2],
                 mass_arrs[0], mass_err_arrs[0], mass_arrs[1], mass_err_arrs[1], mass_arrs[2],
                 mass_err_arrs[2], alpha_B3B6, alpha_B3B6_err, alpha_B6B7, alpha_B6B7_err],
                 names=['D_ID', 'int_flux_B3', 'int_flux_err_B3', 'int_flux_B6', 'int_flux_err_B6',
                         'int_flux_B7', 'int_flux_err_B7', 'inclination_B3', 'inclination_err_B3',
                         'inclination_B6', 'inclination_err_B6', 'inclination_B7', 'inclination_err_B7',
                         'lower_lum_B3', 'lower_lum_B6', 'lower_lum_B7', 'dust_mass_B3', 'dust_mass_err_B3',
                         'dust_mass_B6', 'dust_mass_err_B6', 'dust_mass_B7', 'dust_mass_err_B7', 'alpha_B3B6',
                         'alpha_B3B6_err', 'alpha_B6B7', 'alpha_B6B7_err'],
                 dtype=['i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8',
                        'f8','f8','f8','f8','f8','f8','f8','f8','f8'])    

    calc_tables.append(tab)

data = conv_data
calc_tab = calc_tables[0]

#for src in calc_tables[1]['D_ID']:
#    calc_ind = np.where(calc_tab['D_ID'] == src)[0]
#    calc_tab[calc_ind] = calc_tables[1][np.where(calc_tables[1]['D_ID'] == src)[0]]

#this makes the calculated quantites in the table for sources 3 and 13 from the nonconvolved data

EisnerID = Column(np.array(np.repeat('none', len(data)), dtype='S10'), name='Eisner_ID')

table_meas_misc = Table((data['D_ID'], EisnerID, data['RA_B3'], data['RA_err_B3'], data['DEC_B3'],
                         data['DEC_err_B3'], calc_tab['alpha_B3B6'], calc_tab['alpha_B3B6_err'],
                         calc_tab['alpha_B6B7'], calc_tab['alpha_B6B7_err']))
                         
table_meas_B3 = Table((data['D_ID'], EisnerID, data['ap_flux_B3'], data['ap_flux_err_B3'], data['gauss_amp_B3'],
                       data['gauss_amp_err_B3'], calc_tab['int_flux_B3'],
                       calc_tab['int_flux_err_B3'], data['fwhm_maj_B3'],
                       data['fwhm_maj_err_B3'], data['fwhm_min_B3'],
                       data['fwhm_min_err_B3'], data['pa_B3'], data['pa_err_B3'],
                       data['fwhm_maj_deconv_B3'], data['fwhm_min_deconv_B3'], data['pa_deconv_B3'], 
                       calc_tab['inclination_B3'], calc_tab['inclination_err_B3']))

table_meas_B6 = Table((data['D_ID'], data['ap_flux_B6'],
                       data['ap_flux_err_B6'], data['gauss_amp_B6'], data['gauss_amp_err_B6'],
                       calc_tab['int_flux_B6'], calc_tab['int_flux_err_B6'],
                       data['fwhm_maj_B6'], data['fwhm_maj_err_B6'], data['fwhm_min_B6'],
                       data['fwhm_min_err_B6'], data['pa_B6'], data['pa_err_B6'],
                       data['fwhm_maj_deconv_B6'], data['fwhm_min_deconv_B6'], data['pa_deconv_B6'], 
                       calc_tab['inclination_B6'], calc_tab['inclination_err_B6']))


table_meas_B7 = Table((data['D_ID'], data['ap_flux_B7'], data['ap_flux_err_B7'],
                       data['gauss_amp_B7'], data['gauss_amp_err_B7'],
                       calc_tab['int_flux_B7'], calc_tab['int_flux_err_B7'], data['fwhm_maj_B7'],
                       data['fwhm_maj_err_B7'], data['fwhm_min_B7'],
                       data['fwhm_min_err_B7'], data['pa_B7'], data['pa_err_B7'],
                       data['fwhm_maj_deconv_B7'], data['fwhm_min_deconv_B7'], data['pa_deconv_B7'],
                       calc_tab['inclination_B7'], calc_tab['inclination_err_B7']))


opt_depth = np.repeat('-', len(data['D_ID']))

table_inf = Table((data['D_ID'], calc_tab['lower_lum_B3'], calc_tab['lower_lum_B6'], calc_tab['lower_lum_B7'],
                calc_tab['dust_mass_B3'], calc_tab['dust_mass_err_B3'], calc_tab['dust_mass_B6'],
                calc_tab['dust_mass_err_B6'], calc_tab['dust_mass_B7'], calc_tab['dust_mass_err_B7'],
                opt_depth))


'''
nonconv_srcs = [3,13]
for src in nonconv_srcs:
    ind_tabB6 = np.where(table_meas_B6['D_ID'] == src)[0]
    ind_tabB7 = np.where(table_meas_B7['D_ID'] == src)[0]
    ind_tab2 = np.where(table_inf['D_ID'] == src)[0]
    ind_nonconv = np.where(nonconv_data['D_ID'] == src)[0]
    for col in table_meas_B6.colnames:
        if col in nonconv_data.colnames:
            table_meas_B6[col][ind_tabB6] = nonconv_data[col][ind_nonconv]
    for col in table_meas_B7.colnames:
        if col in nonconv_data.colnames:
            table_meas_B7[col][ind_tabB7] = nonconv_data[col][ind_nonconv]
    for col in table_inf.colnames:
        if col in nonconv_data.colnames:
            table_inf[col][ind_tab2] = nonconv_data[col][ind_nonconv]
'''

eisner_tab = ascii.read('../tables/eisner_tbl.txt', format='tab', delimiter='\t')
eisner_tab.remove_column('remove')

eisner_ra = []
eisner_dec = []
for row in eisner_tab:
    eisner_ra.append(RA_to_deg(float(row['alpha'][0:2]), float(row['alpha'][2:4]), float(row['alpha'][4:])))
    eisner_dec.append(DEC_to_deg(float(row['delta'][0:3]), float(row['delta'][3:5]), float(row['delta'][5:])))

   
coord_tab = SkyCoord(ra=table_meas_misc['RA_B3']*u.degree, dec=table_meas_misc['DEC_B3']*u.degree)

eisner_coord = SkyCoord(ra=eisner_ra*u.degree, dec=eisner_dec*u.degree)

idx, d2d, d3d = eisner_coord.match_to_catalog_sky(coord_tab)
#idx is list of indices of table_meas with locations corresponding to eisner_tab
matches = np.where(d2d.value < 0.5*(1/3600))[0] #matches within 0.1 arcsec
for mat in matches:
    table_meas_misc[idx[mat]]['Eisner_ID'] = eisner_tab[mat]['ID']    

table_meas_B3['Eisner_ID'] = table_meas_misc['Eisner_ID']

eis_coord_tab = Table((eisner_tab['ID'], eisner_ra, eisner_dec),names=('ID', 'RA', 'DEC'))
eis_coord_tab.write('../tables/eis_coord_table.fits', format='fits', overwrite=True)

#tables with table_meas_BX have all sources
B3ind = np.where(np.isnan(table_meas_B3['ap_flux_B3']) == False)[0]
B6ind = np.where(np.isnan(table_meas_B6['ap_flux_B6']) == False)[0]
B7ind = np.where(np.isnan(table_meas_B7['ap_flux_B7']) == False)[0]
allband_ind = np.intersect1d(B6ind, B7ind)

only_B6ind = np.setdiff1d(B6ind, B7ind)
only_B3ind = np.setdiff1d(B3ind, B6ind)

meas_B7 = table_meas_B7
meas_B6 = table_meas_B6
meas_B3 = table_meas_B3
meas_misc = table_meas_misc

meas_misc.write('../tables/table_meas_misc.fits', overwrite=True)
meas_B3.write('../tables/table_meas_B3.fits', overwrite=True)


#now clean up tables to put into latex

radec_err = []
alphaB3B6 = []
alphaB6B7 = []

for row in range(len(meas_misc)):
    ra_err_rnd = round_to_n(meas_misc['RA_err_B3'][row], 1)
    dec_err_rnd = round_to_n(meas_misc['DEC_err_B3'][row], 1)
    radec_err.append('$\pm '+str(ra_err_rnd)+'$, $\pm'+str(dec_err_rnd)+'$')
    if np.isnan(meas_misc['alpha_B3B6'][row]) == True:
        alphaB3B6.append('-')
    else:
        alphaB3B6_rnd, alphaB3B6_err_rnd = rounded(meas_misc['alpha_B3B6'][row], meas_misc['alpha_B3B6_err'][row])
        alphaB3B6.append(str(alphaB3B6_rnd)+'$\pm$'+str(alphaB3B6_err_rnd))
    if np.isnan(meas_misc['alpha_B6B7'][row]) == True:
        alphaB6B7.append('-')
    else:
        alphaB6B7_rnd, alphaB6B7_err_rnd = rounded(meas_misc['alpha_B6B7'][row], meas_misc['alpha_B6B7_err'][row])
        alphaB6B7.append(str(alphaB6B7_rnd)+'$\pm$'+str(alphaB6B7_err_rnd))

coords = SkyCoord(ra=meas_misc['RA_B3']*u.deg, dec=meas_misc['DEC_B3']*u.deg)
coords_str = coords.to_string('hmsdms', sep=':')


tablemisc_latex = Table((meas_misc['D_ID'], meas_misc['Eisner_ID'],
                       coords_str, radec_err,  alphaB3B6, alphaB6B7),
                      names=('Source ID', 'E18 ID', 'Coordinates', '$\\alpha, \\delta$ error',
                             '$\\alpha_{B3\\to B6}$', '$\\alpha_{B6\\to B7}$'))


tablemisc_latex_allband = tablemisc_latex[B7ind]
tablemisc_latex_allband['SED_classification'] = ['-','-','-','-','-','thin','thick','-','thin to thick','-','thin to thick','thick','-','-','thin','thin','-','-']
tablemisc_latex_B3B6 = tablemisc_latex[only_B6ind]
tablemisc_latex_B3B6.remove_column('$\\alpha_{B6\\to B7}$')
tablemisc_latex_B3only = tablemisc_latex[only_B3ind]
tablemisc_latex_B3only.remove_column('$\\alpha_{B3\\to B6}$')
tablemisc_latex_B3only.remove_column('$\\alpha_{B6\\to B7}$')
tablemisc_latex_allband.write('../tables/tablemisc_latex_allband.txt', format='latex', overwrite=True)
tablemisc_latex_B3B6.write('../tables/tablemisc_latex_B3B6.txt', format='latex', overwrite=True)
tablemisc_latex_B3only.write('../tables/tablemisc_latex_B3only.txt', format='latex', overwrite=True)



apfluxB3 = []
intfluxB3 = []
gaussampB3 = []
fwhmmajB3 = []
fwhmmajdeconvB3 = []
fwhmminB3 = []
fwhmmindeconvB3 = []
paB3 = []
padeconvB3 = []
inclB3 =[]


for row in range(len(meas_B3)):
    apflux_rnd, apflux_err_rnd = rounded(meas_B3['ap_flux_B3'][row], meas_B3['ap_flux_err_B3'][row])
    apfluxB3.append(str(apflux_rnd*1000)+'$\pm$'+str(apflux_err_rnd*1000)) #units mJy
    intflux_rnd, intflux_err_rnd = rounded(meas_B3['int_flux_B3'][row], meas_B3['int_flux_err_B3'][row])
    intfluxB3.append(str(intflux_rnd*1000)+'$\pm$'+str(intflux_err_rnd*1000)) #units mJy
    gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B3['gauss_amp_B3'][row], meas_B3['gauss_amp_err_B3'][row])
    gaussampB3.append(str(gaussamp_rnd*1000)+'$\pm$'+str(gaussamp_err_rnd*1000)) #units mJy/beam(?)
    fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B3['fwhm_maj_B3'][row], meas_B3['fwhm_maj_err_B3'][row])
    fwhmmajB3.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
    fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B3['fwhm_min_B3'][row], meas_B3['fwhm_min_err_B3'][row])
    fwhmminB3.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
    pa_rnd, pa_err_rnd = rounded(meas_B3['pa_B3'][row], meas_B3['pa_err_B3'][row])
    paB3.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees
    
    if np.isnan(meas_B3['fwhm_maj_deconv_B3'][row]) == True:
        inclB3.append('-')
        fwhmmajdeconvB3.append('-')
        fwhmmindeconvB3.append('-')
        padeconvB3.append('-')
    else:
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B3['fwhm_maj_deconv_B3'][row], meas_B3['fwhm_maj_err_B3'][row])
        fwhmmajdeconvB3.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B3['fwhm_min_deconv_B3'][row], meas_B3['fwhm_min_err_B3'][row])
        fwhmmindeconvB3.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B3['pa_deconv_B3'][row], meas_B3['pa_err_B3'][row])
        padeconvB3.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB3_rnd, inclB3_err_rnd = rounded(meas_B3['inclination_B3'][row], meas_B3['inclination_err_B3'][row])
        inclB3.append(str(inclB3_rnd)+'$\pm$'+str(inclB3_err_rnd))


tableB3_latex = Table((meas_B3['D_ID'], apfluxB3, intfluxB3, gaussampB3, fwhmmajB3, fwhmminB3, paB3,
                       fwhmmajdeconvB3, fwhmmindeconvB3, padeconvB3, inclB3),
                      names=('Source ID', '$F_{\\text{aperture}, B3}$',
                             '$F_{\\text{integrated}, B3}$','$A_{gaussian, B3}$',
                             '$\\text{FWHM}_{\\text{major}, B3}$', '$\\text{FWHM}_{\\text{minor}, B3}$',
                             '$\\theta_{PA, B3}$', '$\\text{FWHM}_{\\text{major, deconvolved}, B3}$',
                             '$\\text{FWHM}_{\\text{minor, deconvolved}, B3}$', '$\\theta_{PA, \\text{deconvolved}, B3}$',
                             '$\\theta_{\\text{inclination}, B3}$'))


tableB3_latex_allband = tableB3_latex[B7ind]
tableB3_latex_B3B6 = tableB3_latex[only_B6ind]
tableB3_latex_B3only = tableB3_latex[only_B3ind]
tableB3_latex_allband.write('../tables/tableB3_latex_allband.txt', format='latex', overwrite=True)
tableB3_latex_B3B6.write('../tables/tableB3_latex_B3B6.txt', format='latex', overwrite=True)
tableB3_latex_B3only.write('../tables/tableB3_latex_B3only.txt', format='latex', overwrite=True)


apfluxB6 = []
intfluxB6 = []
gaussampB6 = []
fwhmmajB6 = []
fwhmmajdeconvB6 = []
fwhmminB6 = []
fwhmmindeconvB6 = []
paB6 = []
padeconvB6 = []
inclB6 =[]


for row in range(len(meas_B6)):
    if np.isnan(meas_B6['ap_flux_B6'][row]) == False:
        apflux_rnd, apflux_err_rnd = rounded(meas_B6['ap_flux_B6'][row], meas_B6['ap_flux_err_B6'][row])
        apfluxB6.append(str(apflux_rnd*1000)+'$\pm$'+str(apflux_err_rnd*1000)) #units mJy
        intflux_rnd, intflux_err_rnd = rounded(meas_B6['int_flux_B6'][row], meas_B6['int_flux_err_B6'][row])
        intfluxB6.append(str(intflux_rnd*1000)+'$\pm$'+str(intflux_err_rnd*1000)) #units mJy
        gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B6['gauss_amp_B6'][row], meas_B6['gauss_amp_err_B6'][row])
        gaussampB6.append(str(gaussamp_rnd*1000)+'$\pm$'+str(gaussamp_err_rnd*1000)) #units mJy/beam(?)
        fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B6['fwhm_maj_B6'][row], meas_B6['fwhm_maj_err_B6'][row])
        fwhmmajB6.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
        fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B6['fwhm_min_B6'][row], meas_B6['fwhm_min_err_B6'][row])
        fwhmminB6.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
        pa_rnd, pa_err_rnd = rounded(meas_B6['pa_B6'][row], meas_B6['pa_err_B6'][row])
        paB6.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees
    else:
        apfluxB6.append('-')
        intfluxB6.append('-')
        gaussampB6.append('-')
        fwhmmajB6.append('-')
        fwhmminB6.append('-')
        paB6.append('-')
        
    if np.isnan(meas_B6['fwhm_maj_deconv_B6'][row]) == True:
        inclB6.append('-')
        fwhmmajdeconvB6.append('-')
        fwhmmindeconvB6.append('-')
        padeconvB6.append('-')
    else:
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B6['fwhm_maj_deconv_B6'][row], meas_B6['fwhm_maj_err_B6'][row])
        fwhmmajdeconvB6.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B6['fwhm_min_deconv_B6'][row], meas_B6['fwhm_min_err_B6'][row])
        fwhmmindeconvB6.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B6['pa_deconv_B6'][row], meas_B6['pa_err_B6'][row])
        padeconvB6.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB6_rnd, inclB6_err_rnd = rounded(meas_B6['inclination_B6'][row], meas_B6['inclination_err_B6'][row])
        inclB6.append(str(inclB6_rnd)+'$\pm$'+str(inclB6_err_rnd))

tableB6_latex = Table((meas_B6['D_ID'], apfluxB6, intfluxB6, gaussampB6, fwhmmajB6, fwhmminB6, paB6,
                       fwhmmajdeconvB6, fwhmmindeconvB6, padeconvB6, inclB6),
                      names=('Source ID', '$F_{\\text{aperture}, B6}$',
                             '$F_{\\text{integrated}, B6}$','$A_{gaussian, B6}$',
                             '$\\text{FWHM}_{\\text{major}, B6}$', '$\\text{FWHM}_{\\text{minor}, B6}$',
                             '$\\theta_{PA, B6}$', '$\\text{FWHM}_{\\text{major, deconvolved}, B6}$',
                             '$\\text{FWHM}_{\\text{minor, deconvolved}, B6}$', '$\\theta_{PA, \\text{deconvolved}, B6}$',
                             '$\\theta_{\\text{inclination}, B6}$'))


tableB6_latex_allband = tableB6_latex[B7ind]
tableB6_latex_B3B6 = tableB6_latex[only_B6ind]
tableB6_latex_allband.write('../tables/tableB6_latex_allband.txt', format='latex', overwrite=True)
tableB6_latex_B3B6.write('../tables/tableB6_latex_B3B6.txt', format='latex', overwrite=True)



apfluxB7 = []
intfluxB7 = []
gaussampB7 = []
fwhmmajB7 = []
fwhmmajdeconvB7 = []
fwhmminB7 = []
fwhmmindeconvB7 = []
paB7 = []
padeconvB7 = []
inclB7 =[]


for row in range(len(meas_B7)):
    if np.isnan(meas_B7['ap_flux_B7'][row]) == False:
        apflux_rnd, apflux_err_rnd = rounded(meas_B7['ap_flux_B7'][row], meas_B7['ap_flux_err_B7'][row])
        apfluxB7.append(str(apflux_rnd*1000)+'$\pm$'+str(apflux_err_rnd*1000)) #units mJy
        intflux_rnd, intflux_err_rnd = rounded(meas_B7['int_flux_B7'][row], meas_B7['int_flux_err_B7'][row])
        intfluxB7.append(str(intflux_rnd*1000)+'$\pm$'+str(intflux_err_rnd*1000)) #units mJy
        gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B7['gauss_amp_B7'][row], meas_B7['gauss_amp_err_B7'][row])
        gaussampB7.append(str(gaussamp_rnd*1000)+'$\pm$'+str(gaussamp_err_rnd*1000)) #units mJy/beam(?)
        fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B7['fwhm_maj_B7'][row], meas_B7['fwhm_maj_err_B7'][row])
        fwhmmajB7.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
        fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B7['fwhm_min_B7'][row], meas_B7['fwhm_min_err_B7'][row])
        fwhmminB7.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
        pa_rnd, pa_err_rnd = rounded(meas_B7['pa_B7'][row], meas_B7['pa_err_B7'][row])
        paB7.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees
    else:
        apfluxB7.append('-')
        intfluxB7.append('-')
        gaussampB7.append('-')
        fwhmmajB7.append('-')
        fwhmminB7.append('-')
        paB7.append('-')
        
    if np.isnan(meas_B7['fwhm_maj_deconv_B7'][row]) == True:
        inclB7.append('-')
        fwhmmajdeconvB7.append('-')
        fwhmmindeconvB7.append('-')
        padeconvB7.append('-')
    else:
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B7['fwhm_maj_deconv_B7'][row], meas_B7['fwhm_maj_err_B7'][row])
        fwhmmajdeconvB7.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B7['fwhm_min_deconv_B7'][row], meas_B7['fwhm_min_err_B7'][row])
        fwhmmindeconvB7.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B7['pa_deconv_B7'][row], meas_B7['pa_err_B7'][row])
        padeconvB7.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB7_rnd, inclB7_err_rnd = rounded(meas_B7['inclination_B7'][row], meas_B7['inclination_err_B7'][row])
        inclB7.append(str(inclB7_rnd)+'$\pm$'+str(inclB7_err_rnd))

tableB7_latex = Table((meas_B7['D_ID'], apfluxB7, intfluxB7, gaussampB7, fwhmmajB7, fwhmminB7, paB7,
                       fwhmmajdeconvB7, fwhmmindeconvB7, padeconvB7, inclB7),
                      names=('Source ID', '$F_{\\text{aperture}, B7}$',
                             '$F_{\\text{integrated}, B7}$','$A_{gaussian, B7}$',
                             '$\\text{FWHM}_{\\text{major}, B7}$', '$\\text{FWHM}_{\\text{minor}, B7}$',
                             '$\\theta_{PA, B7}$', '$\\text{FWHM}_{\\text{major, deconvolved}, B7}$',
                             '$\\text{FWHM}_{\\text{minor, deconvolved}, B7}$', '$\\theta_{PA, \\text{deconvolved}, B7}$',
                             '$\\theta_{\\text{inclination}, B7}$'))


tableB7_latex_allband = tableB7_latex[B7ind]
tableB7_latex_allband.write('../tables/tableB7_latex_allband.txt', format='latex', overwrite=True)


#ascii.write(table1, '/users/jotter/summer_research_2018/tables/measured_vals_table1.txt', format='latex')
#ascii.write(table2, '/users/jotter/summer_research_2018/tables/inferred_vals_table2.txt', format='latex')


#table_meas_B3.write('/users/jotter/summer_research_2018/tables/measured_vals_B3.fits', overwrite=True)
#table_meas_B6.write('/users/jotter/summer_research_2018/tables/measured_vals_B6.fits', overwrite=True)
#table_meas_B7.write('/users/jotter/summer_research_2018/tables/measured_vals_B7.fits', overwrite=True)


LL_B3 = []
LL_B6 = []
LL_B7 = []
DM_B3 = []
DM_B6 = []
DM_B7 = []

for row in range(len(table_inf)):
        DM_B3_rnd, DM_B3_err_rnd = rounded(table_inf['dust_mass_B3'][row], table_inf['dust_mass_err_B3'][row])
        DM_B3.append(str(DM_B3_rnd)+'$\pm$'+str(DM_B3_err_rnd)) #units mJy
        if np.isnan(table_inf['dust_mass_B6'][row]) == True:
            DM_B6.append('-')
            LL_B6.append('-')
        else:
            DM_B6_rnd, DM_B6_err_rnd = rounded(table_inf['dust_mass_B6'][row], table_inf['dust_mass_err_B6'][row])
            DM_B6.append(str(DM_B6_rnd)+'$\pm$'+str(DM_B6_err_rnd))
            LL_B6.append('>'+str(round_to_n(table_inf['lower_lum_B6'][row],2)))
        if np.isnan(table_inf['dust_mass_B7'][row]) == True:
            DM_B7.append('-')
            LL_B7.append('-')
        else:
            DM_B7_rnd, DM_B7_err_rnd = rounded(table_inf['dust_mass_B7'][row], table_inf['dust_mass_err_B7'][row])
            DM_B7.append(str(DM_B7_rnd)+'$\pm$'+str(DM_B7_err_rnd))
            LL_B7.append('>'+str(round_to_n(table_inf['lower_lum_B7'][row],2)))
        LL_B3.append('>'+str(round_to_n(table_inf['lower_lum_B3'][row],2)))
        
        
tableinf_latex = Table((table_inf['D_ID'], LL_B3, LL_B6, LL_B7, DM_B3, DM_B6, DM_B7),
                      names=('ID','$L_{mm, B3}$','$L_{mm, B6}$','$L_{mm, B7}$','$M_{dust, B3}$',
                             '$M_{dust, B6}$','$M_{dust, B7}$'))


inf_all = tableinf_latex[allband_ind]
inf_B3B6 = tableinf_latex[only_B6ind]
inf_B3B6.remove_columns(['$L_{mm, B7}$','$M_{dust, B7}$'])
inf_B3only = tableinf_latex[only_B3ind]
inf_B3only.remove_columns(['$L_{mm, B6}$','$M_{dust, B6}$','$L_{mm, B7}$','$M_{dust, B7}$'])
inf_all.write('/users/jotter/summer_research_2018/tables/inf_allbands.txt', format='latex', overwrite=True)
inf_B3B6.write('/users/jotter/summer_research_2018/tables/inf_B3B6.txt', format='latex', overwrite=True)
inf_B3only.write('/users/jotter/summer_research_2018/tables/inf_B3only.txt', format='latex', overwrite=True)
