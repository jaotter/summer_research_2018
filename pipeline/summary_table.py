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

data =  Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_nov20_ulim.fits')

calc_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_nov20_calc_vals.fits')

irtab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_nov20_full.fits')


B6nonconv = [16,34,71,80,83]
B6nonconv_ind = [np.where(D_ID_b6==str(b6))[0][0] for b6 in B6nonconv]
B7nonconv = [16,18,33,34,50,71,76,80,81,83]
B7nonconv_ind = [np.where(D_ID_b7==str(b7))[0][0] for b7 in B7nonconv]

D_ID = data['D_ID']

for b6ind in B6nonconv_ind:
    D_ID[b6ind] = f'${D_ID_b6[b6ind]}^*$'
for b7ind in B7nonconv_ind:
    D_ID[b7ind] = f'${D_ID_b7[b7ind]}^\dagger$'
for b6b7 in np.intersect1d(B6nonconv_ind, B7nonconv_ind):
    D_ID[b6b7] = f'${D_ID[b6b7]}^{{*\dagger}}$'

MLLA = np.repeat('-', len(data), dtype='S4')
for mlla_row in irtab:
    data_ind = np.where(data['D_ID'] == mlla_row['D_ID'])[0]
    MLLA[data_ind] = mlla_row['MLLA']+mlla_row[

    
#now clean up tables to put into latex

radec_err = []
alphaB3B6 = []
alphaB6B7 = []

for row in range(len(meas_misc)):
    ra_err_rnd = round_to_n((meas_misc['RA_err_B3'][row]*u.degree).to(u.arcsecond).value, 1)
    dec_err_rnd = round_to_n((meas_misc['DEC_err_B3'][row]*u.degree).to(u.arcsecond).value, 1)
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

coords_str_rnd = []
for crd_str in coords_str:
    split = crd_str.split(':')
    ra_as_str = float(split[2].split(' ')[0])
    other_str = split[2].split(' ')[1]
    ra_as = str(round_to_n(ra_as_str, 5))
    dec_as = str(round_to_n(float(split[-1]), 5))
    while len(ra_as) < 6:
        ra_as = ra_as+'0'
    while len(dec_as) < 6:
        dec_as = dec_as+'0'
    coords_str_rnd.append(f'{split[0]}:{split[1]}:{ra_as} {other_str}:{split[3]}:{dec_as}')

tablemisc_latex = Table((D_ID_b6b7, meas_misc['Eisner_ID'],
                       coords_str_rnd, radec_err,  alphaB3B6, alphaB6B7),
                      names=('Source ID', 'E18 ID', 'Coordinates', '$\\alpha, \\delta$ error',
                             '$\\alpha_{B3\\to B6}$', '$\\alpha_{B6\\to B7}$'))


nonB3only_ind = np.delete(np.arange(len(tablemisc_latex)), only_B3ind)
tablemisc_latex_nonB3only = tablemisc_latex[nonB3only_ind]
tablemisc_latex_allband = tablemisc_latex[B7ind]
#tablemisc_latex_allband['SED_classification'] = ['-','-','-','-','-','thin','thick','-','thin to thick','-','thin to thick','thick','-','-','thin','thin','-','-']
tablemisc_latex_B3B6 = tablemisc_latex[only_B6ind]
tablemisc_latex_B3B6.remove_column('$\\alpha_{B6\\to B7}$')
tablemisc_latex_B3only = tablemisc_latex[only_B3ind]
tablemisc_latex_B3only.remove_column('$\\alpha_{B3\\to B6}$')
tablemisc_latex_B3only.remove_column('$\\alpha_{B6\\to B7}$')
#tablemisc_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_allband.txt', format='latex', overwrite=True)
#tablemisc_latex_B3B6.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_B3B6.txt', format='latex', overwrite=True)
#tablemisc_latex_B3only.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_B3only.txt', format='latex', overwrite=True)
#tablemisc_latex_nonB3only.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_all.txt', format='latex', overwrite=True)
tablemisc_latex.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex.txt', format='latex', overwrite=True)

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
    apflux_rnd, apflux_err_rnd = rounded(meas_B3['ap_flux_B3'][row]*1000, meas_B3['ap_flux_err_B3'][row]*1000)
    apfluxB3.append(str(apflux_rnd)+'$\pm$'+str(apflux_err_rnd)) #units mJy
    intflux_rnd, intflux_err_rnd = rounded(meas_B3['int_flux_B3'][row]*1000, meas_B3['int_flux_err_B3'][row]*1000)
    intfluxB3.append(str(intflux_rnd)+'$\pm$'+str(intflux_err_rnd)) #units mJy
    gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B3['gauss_amp_B3'][row]*1000, meas_B3['gauss_amp_err_B3'][row]*1000)
    gaussampB3.append(str(gaussamp_rnd)+'$\pm$'+str(gaussamp_err_rnd)) #units mJy/beam(?)
    fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B3['fwhm_maj_B3'][row], meas_B3['fwhm_maj_err_B3'][row])
    fwhmmajB3.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
    fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B3['fwhm_min_B3'][row], meas_B3['fwhm_min_err_B3'][row])
    fwhmminB3.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
    pa_rnd, pa_err_rnd = rounded(meas_B3['pa_B3'][row], meas_B3['pa_err_B3'][row])
    paB3.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees
    
    if np.isnan(meas_B3['fwhm_maj_deconv_B3'][row]) == True:
        inclB3.append('-')
        fwhmmindeconvB3.append('-')
        padeconvB3.append('-')
        
        if np.isnan(meas_B3['upper_lim_B3'][row]) == False:
            ulim_as = ((meas_B3['upper_lim_B3'][row]*u.AU / (400*u.pc)).decompose()*u.radian).to(u.arcsecond).value
            ulim_str = f'$<{round_to_n(ulim_as, 1)}$'
            fwhmmajdeconvB3.append(ulim_str)
        else:
            fwhmmajdeconvB3.append('-')
    else:
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B3['fwhm_maj_deconv_B3'][row], meas_B3['fwhm_maj_err_B3'][row])
        fwhmmajdeconvB3.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B3['fwhm_min_deconv_B3'][row], meas_B3['fwhm_min_err_B3'][row])
        fwhmmindeconvB3.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B3['pa_deconv_B3'][row], meas_B3['pa_err_B3'][row])
        padeconvB3.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB3_rnd, inclB3_err_rnd = rounded(meas_B3['inclination_B3'][row], meas_B3['inclination_err_B3'][row])
        inclB3.append(str(inclB3_rnd)+'$\pm$'+str(inclB3_err_rnd))

print(apfluxB3)
        
tableB3_latex = Table((meas_B3['D_ID'], apfluxB3, gaussampB3, fwhmmajB3, fwhmminB3, paB3,
                       fwhmmajdeconvB3, fwhmmindeconvB3, padeconvB3, inclB3),
                      names=('ID', '$F_{\\text{ap}}$','$A_{gaussian}$',
                             '$\\text{FWHM}_{\\text{maj}}$', '$\\text{FWHM}_{\\text{min}}$',
                             '$\\theta$', '$\\text{FWHM}_{\\text{maj, deconv}}$',
                             '$\\text{FWHM}_{\\text{min, deconv}}$', '$\\theta_{\\text{deconv}}$',
                             '$i$'))


tableB3_latex_allband = tableB3_latex[B7ind]
tableB3_latex_B3B6 = tableB3_latex[only_B6ind]
tableB3_latex_B3only = tableB3_latex[only_B3ind]
#tableB3_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_allband.txt', format='latex', overwrite=True)
#tableB3_latex_B3B6.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_B3B6.txt', format='latex', overwrite=True)
#tableB3_latex_B3only.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_B3only.txt', format='latex', overwrite=True)
tableB3_latex.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_all.txt', format='latex', overwrite=True)

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
        apflux_rnd, apflux_err_rnd = rounded(meas_B6['ap_flux_B6'][row]*1000, meas_B6['ap_flux_err_B6'][row]*1000)
        apfluxB6.append(str(apflux_rnd)+'$\pm$'+str(apflux_err_rnd)) #units mJy
        intflux_rnd, intflux_err_rnd = rounded(meas_B6['int_flux_B6'][row]*1000, meas_B6['int_flux_err_B6'][row]*1000)
        intfluxB6.append(str(intflux_rnd)+'$\pm$'+str(intflux_err_rnd)) #units mJy
        gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B6['gauss_amp_B6'][row]*1000, meas_B6['gauss_amp_err_B6'][row]*1000)
        gaussampB6.append(str(gaussamp_rnd)+'$\pm$'+str(gaussamp_err_rnd)) #units mJy/beam(?)
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
        fwhmmindeconvB6.append('-')
        padeconvB6.append('-')
        if np.isnan(meas_B6['upper_lim_B6'][row]) == False:
            ulim_as = ((meas_B6['upper_lim_B6'][row]*u.AU / (400*u.pc)).decompose()*u.radian).to(u.arcsecond).value
            ulim_str = f'$<{round_to_n(ulim_as, 1)}$'
            fwhmmajdeconvB6.append(ulim_str)
        else:
            fwhmmajdeconvB6.append('-')
        
    else:
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B6['fwhm_maj_deconv_B6'][row], meas_B6['fwhm_maj_err_B6'][row])
        fwhmmajdeconvB6.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B6['fwhm_min_deconv_B6'][row], meas_B6['fwhm_min_err_B6'][row])
        fwhmmindeconvB6.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B6['pa_deconv_B6'][row], meas_B6['pa_err_B6'][row])
        padeconvB6.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB6_rnd, inclB6_err_rnd = rounded(meas_B6['inclination_B6'][row], meas_B6['inclination_err_B6'][row])
        inclB6.append(str(inclB6_rnd)+'$\pm$'+str(inclB6_err_rnd))


tableB6_latex = Table((D_ID_b6, apfluxB6, gaussampB6, fwhmmajB6, fwhmminB6, paB6,
                       fwhmmajdeconvB6, fwhmmindeconvB6, padeconvB6, inclB6),
                      names=('ID', '$F_{\\text{ap}}$','$A_{gaussian}$',
                             '$\\text{FWHM}_{\\text{maj}}$', '$\\text{FWHM}_{\\text{min}}$',
                             '$\\theta$', '$\\text{FWHM}_{\\text{maj, deconv}}$',
                             '$\\text{FWHM}_{\\text{min, deconv}}$', '$\\theta_{\\text{deconv}}$',
                             '$i$'))

print(len(meas_B6), len(tableB6_latex))



tableB6_latex_allband = tableB6_latex[B7ind]
tableB6_latex_B3B6 = tableB6_latex[only_B6ind]
#tableB6_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_allband.txt', format='latex', overwrite=True)
#tableB6_latex_B3B6.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_B3B6.txt', format='latex', overwrite=True)
tableB6_latex_all = tableB6_latex[B6ind]
tableB6_latex_all.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_all.txt', format='latex', overwrite=True)

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
        apflux_rnd, apflux_err_rnd = rounded(meas_B7['ap_flux_B7'][row]*1000, meas_B7['ap_flux_err_B7'][row]*1000)
        apfluxB7.append(str(apflux_rnd)+'$\pm$'+str(apflux_err_rnd)) #units mJy
        intflux_rnd, intflux_err_rnd = rounded(meas_B7['int_flux_B7'][row]*1000, meas_B7['int_flux_err_B7'][row]*1000)
        intfluxB7.append(str(intflux_rnd)+'$\pm$'+str(intflux_err_rnd)) #units mJy
        gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B7['gauss_amp_B7'][row]*1000, meas_B7['gauss_amp_err_B7'][row]*1000)
        gaussampB7.append(str(gaussamp_rnd)+'$\pm$'+str(gaussamp_err_rnd)) #units mJy/beam(?)
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
        fwhmmindeconvB7.append('-')
        padeconvB7.append('-')
        
        if np.isnan(meas_B7['upper_lim_B7'][row]) == False:
            ulim_as = ((meas_B7['upper_lim_B7'][row]*u.AU / (400*u.pc)).decompose()*u.radian).to(u.arcsecond).value
            ulim_str = f'$<{round_to_n(ulim_as, 1)}$'
            fwhmmajdeconvB7.append(ulim_str)
        else:
            fwhmmajdeconvB7.append('-')

    else:
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B7['fwhm_maj_deconv_B7'][row], meas_B7['fwhm_maj_err_B7'][row])
        fwhmmajdeconvB7.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B7['fwhm_min_deconv_B7'][row], meas_B7['fwhm_min_err_B7'][row])
        fwhmmindeconvB7.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B7['pa_deconv_B7'][row], meas_B7['pa_err_B7'][row])
        padeconvB7.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB7_rnd, inclB7_err_rnd = rounded(meas_B7['inclination_B7'][row], meas_B7['inclination_err_B7'][row])
        inclB7.append(str(inclB7_rnd)+'$\pm$'+str(inclB7_err_rnd))

tableB7_latex = Table((D_ID_b7, apfluxB7, gaussampB7, fwhmmajB7, fwhmminB7, paB7,
                       fwhmmajdeconvB7, fwhmmindeconvB7, padeconvB7, inclB7),
                      names=('ID', '$F_{\\text{ap}}$','$A_{gaussian}$',
                             '$\\text{FWHM}_{\\text{maj}}$', '$\\text{FWHM}_{\\text{min}}$',
                             '$\\theta$', '$\\text{FWHM}_{\\text{maj, deconv}}$',
                             '$\\text{FWHM}_{\\text{min, deconv}}$', '$\\theta_{\\text{deconv}}$',
                             '$i$'))


tableB7_latex_allband = tableB7_latex[B7ind]
tableB7_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tableB7_latex_allband.txt', format='latex', overwrite=True)


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
        
        
tableinf_latex = Table((D_ID_b6b7, LL_B3, LL_B6, LL_B7, DM_B3, DM_B6, DM_B7),
                      names=('ID','$L_{mm, B3}$','$L_{mm, B6}$','$L_{mm, B7}$','$M_{dust, B3}$',
                             '$M_{dust, B6}$','$M_{dust, B7}$'))


inf_all = tableinf_latex[allband_ind]
inf_B3B6 = tableinf_latex[only_B6ind]
inf_B3B6.remove_columns(['$L_{mm, B7}$','$M_{dust, B7}$'])
inf_B3only = tableinf_latex[only_B3ind]
inf_B3only.remove_columns(['$L_{mm, B6}$','$M_{dust, B6}$','$L_{mm, B7}$','$M_{dust, B7}$'])
#inf_all.write('/home/jotter/nrao/tables/latex_tables/inf_allbands.txt', format='latex', overwrite=True)
#inf_B3B6.write('/home/jotter/nrao/tables/latex_tables/inf_B3B6.txt', format='latex', overwrite=True)
#inf_B3only.write('/home/jotter/nrao/tables/latex_tables/inf_B3only.txt', format='latex', overwrite=True)
tableinf_latex.write('/home/jotter/nrao/tables/latex_tables/inf_all.txt', format='latex', overwrite=True)
