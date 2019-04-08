from astropy.table import Table, hstack
from astropy.io import ascii, fits
import astropy.units as u
import numpy as np
import radio_beam
from astropy import constants


data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_fixed.fits')
bands = ['B3', 'B6', 'B7']
freqs = [98, 223.5, 339.7672758867]*u.GHz
dr = '/lustre/aoc/students/jotter/directory/'
fls = [dr+'OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits', dr+'B6_convolved_r0.5.clean0.05mJy.150mplus.deepmask.image.tt0.pbcor.fits', dr+'B7_convolved_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits']


for b in range(len(bands)):
    fl = fits.open(fls[b])
    imgdata = fl[0].data.squeeze()
    header = fl[0].header
    beam = radio_beam.Beam.from_fits_header(header)
    L_arr = []
    R = (np.sin(beam.minor)*(415*u.pc)).to(u.au)/2
    for row in range(len(data)):
        flux = data['ap_flux_'+bands[b]][row]*u.Jy
        T_B = (flux).to(u.K, beam.jtok_equiv(freqs[b]))
        L = (4 * np.pi * R**2 * constants.sigma_sb * (T_B)**4).to(u.L_sun)
        L_arr.append(L.value)

    data['lum_low_lim_'+bands[b]] = L_arr

data.write('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_apflux_fixed_lumll.fits')

