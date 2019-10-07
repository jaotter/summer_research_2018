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
eis_tab = Table.read('/users/jotter/summer_research_2018/tables/eisner_tbl.txt')
print(eis_tab.info)
print(eis_tab['RA'])
