from astropy.table import Table, hstack
from astropy.io import ascii, fits
import astropy.units as u
import numpy as np
import radio_beam
from astropy import constants

B3data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B3_500klplus_img_catalog_B6_ref.txt')
B6data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B6_500klplus_img_catalog_B6_ref.txt')
B7data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B7_500klplus_img_catalog_B6_ref.txt')

freqs = [98,223.5, 339.7672758867]*u.GHz
catalogs = [B3data, B6data, B7data]
names = ['B3', 'B6', 'B7']

cats = []

for c in range(len(catalogs)):
    fl = fits.open('/lustre/aoc/students/jotter/directory/Orion'+names[c]+'/Orion_SourceI_'+names[c]+'_continuum_r-2.clean0.1mJy.500klplus.deepmask.image.tt0.pbcor.fits')
    imgdata = fl[0].data.squeeze()
    header = fl[0].header
    beam = radio_beam.Beam.from_fits_header(header)
    T_B_arr = []
    L_arr = []
    R = (np.sin(beam.minor)*(415*u.pc)).to(u.au)/2
    for row in range(len(catalogs[c])):
        flux = (catalogs[c]['g_amplitude_r-2.clean0.1mJy.500klplus.deepmask'][row] - catalogs[c]['bg_median_r-2.clean0.1mJy.500klplus.deepmask'][row])*u.Jy
        T_B = (flux).to(u.K, beam.jtok_equiv(freqs[c]))
        L = (4 * np.pi * (R)**2 * constants.sigma_sb * (T_B)**4).to(u.L_sun)
        T_B_arr.append(T_B.value)
        L_arr.append(L.value)
    name = names[c]
    if name == 'B7':
        name = 'B7_hr'
    cat = catalogs[c][['gauss_x_'+name, 'gauss_y_'+name, 'FWHM_major_r-2.clean0.1mJy.500klplus.deepmask', 'major_err_r-2.clean0.1mJy.500klplus.deepmask', 'FWHM_minor_r-2.clean0.1mJy.500klplus.deepmask', 'minor_err_r-2.clean0.1mJy.500klplus.deepmask', 'pa_r-2.clean0.1mJy.500klplus.deepmask', 'pa_err_r-2.clean0.1mJy.500klplus.deepmask', 'g_amplitude_r-2.clean0.1mJy.500klplus.deepmask', 'amp_err_r-2.clean0.1mJy.500klplus.deepmask', 'ap_flux_r-2.clean0.1mJy.500klplus.deepmask', 'ap_flux_err_r-2.clean0.1mJy.500klplus.deepmask', 'bg_median_r-2.clean0.1mJy.500klplus.deepmask', 'bg_ap_r-2.clean0.1mJy.500klplus.deepmask', 'circ_flux_r-2.clean0.1mJy.500klplus.deepmask', 'circ_flux_err_r-2.clean0.1mJy.500klplus.deepmask', 'bg_circ_r-2.clean0.1mJy.500klplus.deepmask']]
    cols = cat.colnames
    for colname in cols:
        col = cat[colname]
        if colname[-8:] == 'deepmask':
            cat.rename_column(colname, colname[:-34]+names[c])

    cat['T_B_'+names[c]] = T_B_arr
    cat['L_'+names[c]] = L_arr
    cats.append(cat)

meta_table = catalogs[0][['D_ID', 'Fb_ID', 'Speak_Fb', 'e_Speak_Fb', 'alpha' , 'e_alpha', 'RA_Fb', 'DEC_Fb', 'e_RAs_Fb', 'e_DEs_Fb', 'HC_ID']]

all_cats = hstack((meta_table, cats[0], cats[1], cats[2]), join_type = 'outer')

all_cats.write('/lustre/aoc/students/jotter/dendro_catalogs/allbands_500klplus.txt', format='ascii')

