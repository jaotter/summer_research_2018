#script to look at sources missed by Eisner 2018
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

data = Table.read('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')

#Eisner FOV limits
eis_RA_upper = 83.836
eis_RA_lower = 83.805
eis_Dec_upper = -5.3715
eis_Dec_lower = -5.402

RA_ind1 = np.where(data['RA_B3'] > eis_RA_lower)[0]
RA_ind2 = np.where(data['RA_B3'] < eis_RA_upper)[0]
RA_ind = np.intersect1d(RA_ind1,RA_ind2)

Dec_ind1 = np.where(data['DEC_B3'] > eis_Dec_lower)[0]
Dec_ind2 = np.where(data['DEC_B3'] < eis_Dec_upper)[0]
Dec_ind = np.intersect1d(Dec_ind1,Dec_ind2)

eis_ind = np.intersect1d(RA_ind, Dec_ind)

#load in eisner data:
eisner_data = Table.read('/users/jotter/summer_research_2018/tables/eisner_tbl.txt', format='ascii')
eis_coord_tab = Table.read('/users/jotter/summer_research_2018/tables/eis_coord_table.fits')

B3_coord = SkyCoord(ra=data['RA_B3']*u.deg, dec=data['DEC_B3']*u.deg)
eis_coord = SkyCoord(ra=eis_coord_tab['RA']*u.deg, dec=eis_coord_tab['DEC']*u.deg)

idx, d2d, d3d = eis_coord.match_to_catalog_sky(B3_coord)
match = np.where(d2d < 0.5*u.arcsec)[0] #match within 0.1 arcsec

data['Eis_ID'] = np.repeat(np.array('-',dtype='S8'), len(data))

for eis_ind in match:
    data['Eis_ID'][idx[eis_ind]] = eis_coord_tab['ID'][eis_ind]

eis_match_ind = np.where(data['Eis_ID'] != '-')[0]
print(len(eis_match_ind))

fov_src_hist, bins = np.histogram(data['fwhm_maj_B3'][eis_ind], density=False)
matched_hist, b = np.histogram(data['fwhm_maj_B3'][eis_match_ind], bins, density=False)

plotpts = []
widths = []
for b in range(len(bins[:-1])): #creating points to plot - midpoints of bins
    plotpts.append(bins[b] + (bins[b+1]-bins[b])/2)
    widths.append((bins[b+1]-bins[b]))
                            
