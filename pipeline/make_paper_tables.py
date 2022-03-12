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

data =  Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')

calc_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_may21_calc_vals_mask_alpha_ulim.fits')

irtab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')

MLLA = Column(np.array(np.repeat('-', len(data)), dtype='S10'), name='MLLA')

match = Table.read('/home/jotter/nrao/summer_research_2018/tables/COUP_Forbrich16_r0.5_may21_full.fits')

for ir_row in irtab:
   
    if ir_row['m_MLLA'] == ' ':
        mlla = str(ir_row['MLLA'])
    else:
        mlla = str(ir_row['MLLA']) + str(ir_row['m_MLLA']).lower()
    ind = np.where(data['ID'] == ir_row['ID'])[0]
    MLLA[ind] = mlla

table_meas_misc = Table((data['ID'], MLLA, match['Seq'], match['COUP_1'], data['RA_B3'], data['RA_err_B3'], data['DEC_B3'],
                         data['DEC_err_B3'], calc_tab['alpha_B3B6'], calc_tab['alpha_B3B6_err'],
                         calc_tab['alpha_B6B7'], calc_tab['alpha_B6B7_err'], calc_tab['alpha_fit'], calc_tab['alpha_fit_err'], calc_tab['alpha_ulim_B3B6'],
                         calc_tab['alpha_B3eis'], calc_tab['alpha_B3eis_err']))
                         
table_meas_B3 = Table((data['ID'], data['ap_flux_B3'], data['ap_flux_err_B3'], data['gauss_amp_B3'],
                       data['gauss_amp_err_B3'], calc_tab['int_flux_B3'],
                       calc_tab['int_flux_err_B3'], data['fwhm_maj_B3'],
                       data['fwhm_maj_err_B3'], data['fwhm_min_B3'],
                       data['fwhm_min_err_B3'], data['pa_B3'], data['pa_err_B3'],
                       data['fwhm_maj_deconv_B3'], data['fwhm_min_deconv_B3'], data['pa_deconv_B3'], 
                       calc_tab['inclination_B3'], calc_tab['inclination_err_B3'], data['upper_lim_B3'], calc_tab['dust_mass_B3']))

table_meas_B6 = Table((data['ID'], data['ap_flux_B6'],
                       data['ap_flux_err_B6'], data['gauss_amp_B6'], data['gauss_amp_err_B6'],
                       calc_tab['int_flux_B6'], calc_tab['int_flux_err_B6'],
                       data['fwhm_maj_B6'], data['fwhm_maj_err_B6'], data['fwhm_min_B6'],
                       data['fwhm_min_err_B6'], data['pa_B6'], data['pa_err_B6'],
                       data['fwhm_maj_deconv_B6'], data['fwhm_min_deconv_B6'], data['pa_deconv_B6'], 
                       calc_tab['inclination_B6'], calc_tab['inclination_err_B6'], data['upper_lim_B6'], calc_tab['dust_mass_B6'],
                       data['B6_flux_ulim'], calc_tab['dust_mass_ulim_B6']))


table_meas_B7 = Table((data['ID'], data['ap_flux_B7'], data['ap_flux_err_B7'],
                       data['gauss_amp_B7'], data['gauss_amp_err_B7'],
                       calc_tab['int_flux_B7'], calc_tab['int_flux_err_B7'], data['fwhm_maj_B7'],
                       data['fwhm_maj_err_B7'], data['fwhm_min_B7'],
                       data['fwhm_min_err_B7'], data['pa_B7'], data['pa_err_B7'],
                       data['fwhm_maj_deconv_B7'], data['fwhm_min_deconv_B7'], data['pa_deconv_B7'],
                       calc_tab['inclination_B7'], calc_tab['inclination_err_B7'], data['upper_lim_B7'], calc_tab['dust_mass_B7'],
                       data['B7_flux_ulim'], calc_tab['dust_mass_ulim_B7']))



table_inf = Table((data['ID'], calc_tab['lower_lum_B3'], calc_tab['lower_lum_B6'], calc_tab['lower_lum_B7'],
                calc_tab['dust_mass_B3'], calc_tab['dust_mass_err_B3'], calc_tab['dust_mass_B6'],
                calc_tab['dust_mass_err_B6'], calc_tab['dust_mass_B7'], calc_tab['dust_mass_err_B7']))


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

#meas_misc.write('/home/jotter/nrao/tables/table_meas_misc.fits', overwrite=True)
#meas_B3.write('/home/jotter/nrao/tables/table_meas_B3.fits', overwrite=True)

D_ID = np.array(meas_misc['ID'], dtype="<U20")
D_ID_b6b7 = np.array(meas_misc['ID'], dtype="<U20") #include marks for b6 and b7 nonconv srcs
D_ID_b6 = np.array(meas_misc['ID'], dtype="<U20") #only include b6 nonconv srcs
D_ID_b7 = np.array(meas_misc['ID'], dtype="<U20") #only include b7 nonconv srcs

#B6nonconv = [14,30,62,67,69] #updated - these are Seq
B6nonconv = [30,36,45,56,63,68,70,71] #ID values, updated may21
B6nonconv_ind = [np.where(D_ID_b6==str(b6))[0][0] for b6 in B6nonconv]
#B7nonconv = [14,16,29,30,45,62,65,67,68,69] # updated to Seq [16,18,33,34,50,71,76,80,81,83]
B7nonconv = [14,16,29,30,45,63,66,68,69,70] #ID values 
B7nonconv_ind = [np.where(D_ID_b7==str(b7))[0][0] for b7 in B7nonconv]

for b6ind in B6nonconv_ind:
    D_ID_b6[b6ind] = f'{D_ID_b6[b6ind]}$^b$'
    D_ID_b6b7[b6ind] = f'{D_ID[b6ind]}$^b$'
for b7ind in B7nonconv_ind:
    D_ID_b7[b7ind] = f'{D_ID_b7[b7ind]}$^c$'
    D_ID_b6b7[b7ind] = f'{D_ID[b7ind]}$^c$'
for b6b7 in np.intersect1d(B6nonconv_ind, B7nonconv_ind):
    D_ID_b6b7[b6b7] = f'{D_ID[b6b7]}$^{{bc}}$' 

D_ID_full = D_ID_b6b7
#new_srcs = [8,10,32,33,50,53,63,70,74,75,79,115]#updated to new sources [9,10,22,24,36,37,39,45,55,59,61,71,72,79,81,84]
new_srcs = [8, 10, 32, 33, 50, 54, 64, 71, 75, 76, 80, 118, 119, 123, 124]
for new in new_srcs:
    new_ind = np.where(D_ID == str(new))[0]
    D_ID_full[new_ind] = f'{D_ID_full[new_ind[0]]}$^a$'


eisb7_ind = np.where(np.isnan(calc_tab['alpha_B3eis']) == False)[0]
for eis in eisb7_ind:
    D_ID_full[eis] = f'{D_ID_full[eis]}$^d$'

#now clean up tables to put into latex

radec_err = []
alphaB3B6 = []
alphaB6B7 = []
alphafit = []

for row in range(len(meas_misc)):
    ra_err_rnd = round_to_n((meas_misc['RA_err_B3'][row]*u.degree).to(u.arcsecond).value, 1)
    dec_err_rnd = round_to_n((meas_misc['DEC_err_B3'][row]*u.degree).to(u.arcsecond).value, 1)
    radec_err.append('$\pm '+str(ra_err_rnd)+'$, $\pm'+str(dec_err_rnd)+'$')
    if np.isnan(meas_misc['alpha_B3B6'][row]) == True:
        if np.isnan(meas_misc['alpha_ulim_B3B6'][row]) == True:
            alphaB3B6.append('-')
        else:
            alphaB3B6_ulim_rnd, err_rnd = rounded(meas_misc['alpha_ulim_B3B6'][row], np.abs(meas_misc['alpha_ulim_B3B6'][row]/10), extra=0)
            alphaB3B6.append('$<$'+str(alphaB3B6_ulim_rnd))
    else:
        alphaB3B6_rnd, alphaB3B6_err_rnd = rounded(meas_misc['alpha_B3B6'][row], meas_misc['alpha_B3B6_err'][row], extra=0)
        alphaB3B6.append(str(alphaB3B6_rnd)+'$\pm$'+str(alphaB3B6_err_rnd))
    if np.isnan(meas_misc['alpha_B6B7'][row]) == True:
        alphaB6B7.append('-')
    else:
        alphaB6B7_rnd, alphaB6B7_err_rnd = rounded(meas_misc['alpha_B6B7'][row], meas_misc['alpha_B6B7_err'][row], extra=0)
        alphaB6B7.append(str(alphaB6B7_rnd)+'$\pm$'+str(alphaB6B7_err_rnd))
    if np.isnan(meas_misc['alpha_fit'][row]) == True:
        if np.isnan(meas_misc['alpha_B3eis'][row]) == False:
            b3eis_rnd, b3eis_err_rnd = rounded(meas_misc['alpha_B3eis'][row], meas_misc['alpha_B3eis_err'][row], extra=0)
            alphafit.append(str(b3eis_rnd)+'$\pm$'+str(b3eis_err_rnd))
        elif alphaB3B6[-1] != '-':
            alphafit.append(alphaB3B6[-1])
        elif np.isnan(meas_misc['alpha_B3eis'][row]) == True and alphaB3B6[-1] == '-':
            alphafit.append('-')
    else:
        alphafit_rnd, alphafit_err_rnd = rounded(meas_misc['alpha_fit'][row], meas_misc['alpha_fit_err'][row], extra=0)
        alphafit.append(str(alphafit_rnd)+'$\pm$'+str(alphafit_err_rnd))
                              
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

#make this table cover all sources
tablemisc_latex = Table((D_ID_full, meas_misc['MLLA'], meas_misc['Seq'], meas_misc['COUP_1'],
                         coords_str_rnd, radec_err,  alphaB3B6, alphaB6B7, alphafit),
                      names=('ID', 'MLLA', 'Forbrich2016', 'COUP', 'Coordinates', '$\\alpha, \\delta$ error',
                             '$\\alpha_{B3\\to B6}$', '$\\alpha_{B6\\to B7}$', '$\\alpha^{(1)}$'))


nonB3only_ind = np.delete(np.arange(len(tablemisc_latex)), only_B3ind)
tablemisc_latex_nonB3only = tablemisc_latex[nonB3only_ind]
tablemisc_latex_allband = tablemisc_latex[B7ind]
tablemisc_latex_B3B6 = tablemisc_latex[only_B6ind]
tablemisc_latex_B3B6.remove_column('$\\alpha_{B6\\to B7}$')
tablemisc_latex_B3only = tablemisc_latex[only_B3ind]
tablemisc_latex_B3only.remove_column('$\\alpha_{B3\\to B6}$')
tablemisc_latex_B3only.remove_column('$\\alpha_{B6\\to B7}$')
#tablemisc_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_allband.txt', format='latex', overwrite=True)
#tablemisc_latex_B3B6.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_B3B6.txt', format='latex', overwrite=True)
#tablemisc_latex_B3only.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_B3only.txt', format='latex', overwrite=True)
#tablemisc_latex_nonB3only.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_all.txt', format='latex', overwrite=True)

print(tablemisc_latex['ID'][0])
tablemisc_latex['ID'][0] = '127$^d$'
tablemisc_latex.add_row(tablemisc_latex[0])
tablemisc_latex.remove_row(0)

tablemisc_latex.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex.txt', format='latex', overwrite=True)
tablemisc_latex.write('/home/jotter/nrao/tables/latex_tables/final_tables/table4_full.fits', format='fits', overwrite=True)

trunc_length = 6
tablemisc_latex_trunc = tablemisc_latex[:trunc_length]
tablemisc_latex_trunc.write('/home/jotter/nrao/tables/latex_tables/tablemisc_latex_trunc.txt', format='latex', overwrite=True)


apfluxB3 = []
gaussampB3 = []
fwhmmajB3 = []
fwhmmajdeconvB3 = []
fwhmminB3 = []
fwhmmindeconvB3 = []
paB3 = []
padeconvB3 = []
inclB3 =[]
DM_B3 = []

for row in range(len(meas_B3)):
    apflux_rnd, apflux_err_rnd = rounded(meas_B3['ap_flux_B3'][row]*1000, meas_B3['ap_flux_err_B3'][row]*1000, extra=0)
    apfluxB3.append(str(apflux_rnd)+'$\pm$'+str(apflux_err_rnd)) #units mJy
    gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B3['gauss_amp_B3'][row]*1000, meas_B3['gauss_amp_err_B3'][row]*1000, extra=0)
    gaussampB3.append(str(gaussamp_rnd)+'$\pm$'+str(gaussamp_err_rnd)) #units mJy/beam(?)
    fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B3['fwhm_maj_B3'][row], meas_B3['fwhm_maj_err_B3'][row], extra=0)
    fwhmmajB3.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
    fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B3['fwhm_min_B3'][row], meas_B3['fwhm_min_err_B3'][row], extra=0)
    fwhmminB3.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
    pa_rnd, pa_err_rnd = rounded(meas_B3['pa_B3'][row], meas_B3['pa_err_B3'][row], extra=0)
    paB3.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees
    DM_B3_rnd, DM_B3_err_rnd = rounded(table_inf['dust_mass_B3'][row], table_inf['dust_mass_err_B3'][row], extra=0)
    DM_B3.append(str(strip_trailing_zeros(str(DM_B3_rnd)))+'$\pm$'+str(strip_trailing_zeros(str(DM_B3_err_rnd)))) #units mJy

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
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B3['fwhm_maj_deconv_B3'][row], meas_B3['fwhm_maj_err_B3'][row], extra=0)
        fwhmmajdeconvB3.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B3['fwhm_min_deconv_B3'][row], meas_B3['fwhm_min_err_B3'][row], extra=0)
        fwhmmindeconvB3.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B3['pa_deconv_B3'][row], meas_B3['pa_err_B3'][row], extra=0)
        padeconvB3.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB3_rnd, inclB3_err_rnd = rounded(meas_B3['inclination_B3'][row], meas_B3['inclination_err_B3'][row], extra=0)
        inclB3.append(strip_trailing_zeros(str(inclB3_rnd))+'$\pm$'+strip_trailing_zeros(str(inclB3_err_rnd)))


poorfit_B3 = [113,114,124,125]
B3_ID = [f'{b3id}' for b3id in meas_B3['ID']]
for pf in poorfit_B3:
   B3_ID[pf] = B3_ID[pf]+'$^*$'

B3_ID = np.array(B3_ID)


tableB3_latex = Table((B3_ID, apfluxB3, gaussampB3, fwhmmajB3, fwhmminB3, paB3,
                       fwhmmajdeconvB3, fwhmmindeconvB3, padeconvB3, inclB3, DM_B3),
                      names=('ID', '$F_{\\text{ap}}$','$A_{gaussian}$',
                             '$\\text{FWHM}_{\\text{maj}}$', '$\\text{FWHM}_{\\text{min}}$',
                             '$\\theta$', '$\\text{FWHM}_{\\text{maj, deconv}}$',
                             '$\\text{FWHM}_{\\text{min, deconv}}$', '$\\theta_{\\text{deconv}}$',
                             '$i$', '$M_{dust, B3}$'))


print(tableB3_latex['ID'][0])
tableB3_latex['ID'][0] = 127
tableB3_latex.add_row(tableB3_latex[0])
tableB3_latex.remove_row(0)

#tableB3_latex_allband = tableB3_latex[B7ind]
#tableB3_latex_B3B6 = tableB3_latex[only_B6ind]
#tableB3_latex_B3only = tableB3_latex[only_B3ind]
#tableB3_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_allband.txt', format='latex', overwrite=True)
#tableB3_latex_B3B6.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_B3B6.txt', format='latex', overwrite=True)
#tableB3_latex_B3only.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_B3only.txt', format='latex', overwrite=True)
tableB3_latex.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_all.txt', format='latex', overwrite=True)
tableB3_latex.write('/home/jotter/nrao/tables/latex_tables/final_tables/table5_full.fits', format='fits', overwrite=True)


trunc_length = 6
tableB3_latex_trunc = tableB3_latex[:trunc_length]
tableB3_latex_trunc.write('/home/jotter/nrao/tables/latex_tables/tableB3_latex_all_trunc.txt', format='latex', overwrite=True)


apfluxB6 = []
gaussampB6 = []
fwhmmajB6 = []
fwhmmajdeconvB6 = []
fwhmminB6 = []
fwhmmindeconvB6 = []
paB6 = []
padeconvB6 = []
inclB6 =[]
DM_B6 = []


for row in range(len(meas_B3)):
    if np.isnan(meas_B6['fwhm_maj_B6'][row]) == False:
        gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B6['gauss_amp_B6'][row]*1000, meas_B6['gauss_amp_err_B6'][row]*1000, extra=0)
        gaussampB6.append(str(gaussamp_rnd)+'$\pm$'+str(gaussamp_err_rnd)) #units mJy/beam(?)
        fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B6['fwhm_maj_B6'][row], meas_B6['fwhm_maj_err_B6'][row], extra=0)
        fwhmmajB6.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
        fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B6['fwhm_min_B6'][row], meas_B6['fwhm_min_err_B6'][row], extra=0)
        fwhmminB6.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
        pa_rnd, pa_err_rnd = rounded(meas_B6['pa_B6'][row], meas_B6['pa_err_B6'][row], extra=0)
        paB6.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees

        apflux_rnd, apflux_err_rnd = rounded(meas_B6['ap_flux_B6'][row]*1000, meas_B6['ap_flux_err_B6'][row]*1000, extra=0)
        apfluxB6.append(str(apflux_rnd)+'$\pm$'+str(apflux_err_rnd)) #units mJy
        DM_B6_rnd, DM_B6_err_rnd = rounded(table_inf['dust_mass_B6'][row], table_inf['dust_mass_err_B6'][row], extra=0)
        DM_B6.append(strip_trailing_zeros(str(DM_B6_rnd))+'$\pm$'+strip_trailing_zeros(str(DM_B6_err_rnd)))
        
    else:
        gaussampB6.append('-')
        fwhmmajB6.append('-')
        fwhmminB6.append('-')
        paB6.append('-')

        if np.isnan(meas_B6['B6_flux_ulim'][row]) == False:
            apflux_ulim_rnd, small = rounded(meas_B6['B6_flux_ulim'][row]*1000,meas_B6['B6_flux_ulim'][row]*10, extra=0)
            apfluxB6.append(f'$<${apflux_ulim_rnd}')
            DM_ulim_rnd, small = rounded(meas_B6['dust_mass_ulim_B6'][row], meas_B6['dust_mass_ulim_B6'][row]/100, extra=0)
            DM_B6.append(f'$<${DM_ulim_rnd}')

            
        else:
            apfluxB6.append('-')
            DM_B6.append('-')

            
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
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B6['fwhm_maj_deconv_B6'][row], meas_B6['fwhm_maj_err_B6'][row], extra=0)
        fwhmmajdeconvB6.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B6['fwhm_min_deconv_B6'][row], meas_B6['fwhm_min_err_B6'][row], extra=0)
        fwhmmindeconvB6.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B6['pa_deconv_B6'][row], meas_B6['pa_err_B6'][row], extra=0)
        padeconvB6.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB6_rnd, inclB6_err_rnd = rounded(meas_B6['inclination_B6'][row], meas_B6['inclination_err_B6'][row], extra=0)
        inclB6.append(strip_trailing_zeros(str(inclB6_rnd))+'$\pm$'+strip_trailing_zeros(str(inclB6_err_rnd)))


tableB6_latex = Table((D_ID_b6, apfluxB6, gaussampB6, fwhmmajB6, fwhmminB6, paB6,
                       fwhmmajdeconvB6, fwhmmindeconvB6, padeconvB6, inclB6, DM_B6),
                      names=('ID', '$F_{\\text{ap}}$','$A_{gaussian}$',
                             '$\\text{FWHM}_{\\text{maj}}$', '$\\text{FWHM}_{\\text{min}}$',
                             '$\\theta$', '$\\text{FWHM}_{\\text{maj, deconv}}$',
                             '$\\text{FWHM}_{\\text{min, deconv}}$', '$\\theta_{\\text{deconv}}$',
                             '$i$', '$M_{dust, B6}$'))


#tableB6_latex_allband = tableB6_latex[B7ind]
#tableB6_latex_B3B6 = tableB6_latex[only_B6ind]
#tableB6_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_allband.txt', format='latex', overwrite=True)
#tableB6_latex_B3B6.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_B3B6.txt', format='latex', overwrite=True)
ind_b6 = np.where(tableB6_latex['$F_{\\text{ap}}$'] != '-')
tableB6_latex_all = tableB6_latex[ind_b6]
tableB6_latex_all.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_all.txt', format='latex', overwrite=True)
tableB6_latex_all.write('/home/jotter/nrao/tables/latex_tables/final_tables/table6_full.fits', format='fits', overwrite=True)


trunc_length = 6
tableB6_latex_all_trunc = tableB6_latex_all[:trunc_length]
tableB6_latex_all_trunc.write('/home/jotter/nrao/tables/latex_tables/tableB6_latex_all_trunc.txt', format='latex', overwrite=True)

apfluxB7 = []
gaussampB7 = []
fwhmmajB7 = []
fwhmmajdeconvB7 = []
fwhmminB7 = []
fwhmmindeconvB7 = []
paB7 = []
padeconvB7 = []
inclB7 =[]
DM_B7 = []


for row in range(len(meas_B7)):
    if np.isnan(meas_B7['ap_flux_B7'][row]) == False:
        gaussamp_rnd, gaussamp_err_rnd = rounded(meas_B7['gauss_amp_B7'][row]*1000, meas_B7['gauss_amp_err_B7'][row]*1000, extra=0)
        gaussampB7.append(str(gaussamp_rnd)+'$\pm$'+str(gaussamp_err_rnd)) #units mJy/beam(?)
        fwhmmaj_rnd, fwhmmaj_err_rnd = rounded(meas_B7['fwhm_maj_B7'][row], meas_B7['fwhm_maj_err_B7'][row], extra=0)
        fwhmmajB7.append(str(fwhmmaj_rnd)+'$\pm$'+str(fwhmmaj_err_rnd)) #units arcseconds
        fwhmmin_rnd, fwhmmin_err_rnd = rounded(meas_B7['fwhm_min_B7'][row], meas_B7['fwhm_min_err_B7'][row], extra=0)
        fwhmminB7.append(str(fwhmmin_rnd)+'$\pm$'+str(fwhmmin_err_rnd)) #units arcseconds
        pa_rnd, pa_err_rnd = rounded(meas_B7['pa_B7'][row], meas_B7['pa_err_B7'][row], extra=0)
        paB7.append(str(pa_rnd)+'$\pm$'+str(pa_err_rnd)) #units degrees
        apflux_rnd, apflux_err_rnd = rounded(meas_B7['ap_flux_B7'][row]*1000, meas_B7['ap_flux_err_B7'][row]*1000, extra=0)
        apfluxB7.append(str(apflux_rnd)+'$\pm$'+str(apflux_err_rnd)) #units mJy
        DM_B7_rnd, DM_B7_err_rnd = rounded(table_inf['dust_mass_B7'][row], table_inf['dust_mass_err_B7'][row], extra=0)
        DM_B7.append(strip_trailing_zeros(str(DM_B7_rnd))+'$\pm$'+strip_trailing_zeros(str(DM_B7_err_rnd)))

    else:
        gaussampB7.append('-')
        fwhmmajB7.append('-')
        fwhmminB7.append('-')
        paB7.append('-')
        
        if np.isnan(meas_B7['B7_flux_ulim'][row]) == False:
            apflux_ulim_rnd, small = rounded(meas_B7['B7_flux_ulim'][row]*1000,meas_B7['B7_flux_ulim'][row]*10, extra=0)
            apfluxB7.append(f'$<${apflux_ulim_rnd}')
            DM_ulim_rnd, small = rounded(meas_B7['dust_mass_ulim_B7'][row], meas_B7['dust_mass_ulim_B7'][row]/100, extra=0)
            DM_B7.append(f'$<${DM_ulim_rnd}')

        else:
            apfluxB7.append('-')
            DM_B7.append('-')

            
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
        fwhmmajd_rnd, fwhmmajd_err_rnd = rounded(meas_B7['fwhm_maj_deconv_B7'][row], meas_B7['fwhm_maj_err_B7'][row], extra=0)
        fwhmmajdeconvB7.append(str(fwhmmajd_rnd)+'$\pm$'+str(fwhmmajd_err_rnd)) #units arcseconds
        fwhmmind_rnd, fwhmmind_err_rnd = rounded(meas_B7['fwhm_min_deconv_B7'][row], meas_B7['fwhm_min_err_B7'][row], extra=0)
        fwhmmindeconvB7.append(str(fwhmmind_rnd)+'$\pm$'+str(fwhmmind_err_rnd)) #units arcseconds
        pad_rnd, pad_err_rnd = rounded(meas_B7['pa_deconv_B7'][row], meas_B7['pa_err_B7'][row], extra=0)
        padeconvB7.append(str(pad_rnd)+'$\pm$'+str(pad_err_rnd)) #units degrees
        inclB7_rnd, inclB7_err_rnd = rounded(meas_B7['inclination_B7'][row], meas_B7['inclination_err_B7'][row], extra=0)
        inclB7.append(strip_trailing_zeros(str(inclB7_rnd))+'$\pm$'+strip_trailing_zeros(str(inclB7_err_rnd)))

tableB7_latex = Table((D_ID_b7, apfluxB7, gaussampB7, fwhmmajB7, fwhmminB7, paB7,
                       fwhmmajdeconvB7, fwhmmindeconvB7, padeconvB7, inclB7, DM_B7),
                      names=('ID', '$F_{\\text{ap}}$','$A_{gaussian}$',
                             '$\\text{FWHM}_{\\text{maj}}$', '$\\text{FWHM}_{\\text{min}}$',
                             '$\\theta$', '$\\text{FWHM}_{\\text{maj, deconv}}$',
                             '$\\text{FWHM}_{\\text{min, deconv}}$', '$\\theta_{\\text{deconv}}$',
                             '$i$','$M_{dust, B7}$'))


ind_b7 = np.where(tableB7_latex['$F_{\\text{ap}}$'] != '-')
tableB7_latex_allband = tableB7_latex[ind_b7]
tableB7_latex_allband.write('/home/jotter/nrao/tables/latex_tables/tableB7_latex_allband.txt', format='latex', overwrite=True)
tableB7_latex_allband.write('/home/jotter/nrao/tables/latex_tables/final_tables/table7_full.fits', format='fits', overwrite=True)


trunc_length = 6
tableB7_latex_allband_trunc = tableB7_latex_allband[:trunc_length]
tableB7_latex_allband_trunc.write('/home/jotter/nrao/tables/latex_tables/tableB7_latex_allband_trunc.txt', format='latex', overwrite=True)
