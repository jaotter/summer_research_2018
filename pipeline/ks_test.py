from astropy.table import Table
from astropy.io import fits, ascii
import numpy as np
from scipy import stats
import astropy.units as u
import radio_beam


B3img = fits.open('/lustre/aoc/students/jotter/directory/OrionB3/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
B3head = B3img[0].header
B3beam = radio_beam.Beam.from_fits_header(B3head)
B3beamau = ((B3beam.major.to(u.rad))*(414*u.pc).to(u.au)).value

data = Table.read('../tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')
eisner = ascii.read('../tables/eisner_tbl.txt',format='tab')

eisner_R = eisner['R_disk']
data_R = data['fwhm_maj_deconv_B3']*u.arcsec

eis_meas_ind = np.where(eisner['R_disk'] != '<5')[0]
eisner_Redit = eisner['R_disk']
eis_UL_ind = np.where(eisner['R_disk'] == '<5')[0]
eisner_Redit[eis_UL_ind] = '5.0'
eisner_Redit = [float(x.split()[0]) for x in eisner['R_disk']]

eisner_R = [float(x.split()[0]) for x in eisner['R_disk'][eis_meas_ind]]


deconv_ind = np.where(np.isnan(data_R) == False)[0]
data_R = (data_R[deconv_ind].to(u.rad)*(414*u.pc).to(u.au)).value

print(B3beamau)
print(stats.ks_2samp(eisner_R, data_R))
