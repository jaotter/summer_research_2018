#SIMBAD name [OW94] 163-224
#this match might not be real - on SIMBAD HC2000 data and OW94 not officially matched

import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'))

OW94 = ascii.read('/lustre/aoc/students/jotter/tables/OW94.txt', data_start=7, data_end=387, header_start=4)
OW94['RA_OW94'] = RA_to_deg(OW94['RAh'], OW94['RAm'], OW94['RAs'])
OW94['DEC_OW94'] = DEC_to_deg(OW94['DEd'], OW94['DEm'], OW94['DEs'])
OW94.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
OW94['D_ID'] = np.full(len(OW94), -1)

OW_coord = SkyCoord(OW94['RA_OW94']*u.degree, OW94['DEC_OW94']*u.degree)
B3_coord = SkyCoord(data['gauss_x_B3']*u.degree, data['gauss_y_B3']*u.degree)

idx, d2d, d3d = OW_coord.match_to_catalog_sky(B3_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	OW94[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

OW_joined = join(data, OW94, keys='D_ID', join_type='left')

ind83 = np.where(OW_joined['_idx_B3'] == 83)
ind_ref = np.where(OW_joined['_idx_B3'] == 82)

HC_ra_err = (0.017*u.arcsecond).to(u.degree)
HC_dec_err = (0.028*u.arcsecond).to(u.degree)

RAs = [float(OW_joined['RA_OW94'][ind83] - OW_joined['RA_OW94'][ind_ref]),float(OW_joined['RA_HC'][ind83] - OW_joined['RA_HC'][ind_ref]),float(OW_joined['gauss_x_B3'][ind83] - OW_joined['gauss_x_B3'][ind_ref])]
DECs = [float(OW_joined['DEC_OW94'][ind83] - OW_joined['DEC_OW94'][ind_ref]),float(OW_joined['DEC_HC'][ind83] - OW_joined['DEC_HC'][ind_ref]),float(OW_joined['gauss_y_B3'][ind83] - OW_joined['gauss_y_B3'][ind_ref])]
RA_err = [0.1/(3600),HC_ra_err.value,OW_joined['x_err_B3'][ind83]]
DEC_err = [0.1/(3600),HC_dec_err.value,OW_joined['y_err_B3'][ind83]]

times=[0,5+9/12,23+10/12]
times_err=[0.01,0.002, 0.001]


OW94_date = (1,1994)
OW94_err = 0.5/12
HC_date = (10,1999)
HC_err = 0.5/12
B3_date = (10,2017)
B3_err = 0.5/12

dates = [OW94_date,HC_date,B3_date]
times = calc_times(dates)

times_err = [OW94_err,HC_err,B3_err]


v = calc_pm(RAs, DECs, RA_err, DEC_err, times, times_err)
plt.show()
'''
fig = plt.figure()
ax = plt.axes()
#plot trajectory of source 37
OW94 = [5, 35, 16.29, -5, 22,24.1]
OW94_RA = RA_to_deg(OW94[0], OW94[1], OW94[2])
OW94_DEC = DEC_to_deg(OW94[3], OW94[4], OW94[5])
print(OW94_RA, OW94_DEC)
plt.scatter([OW94_RA, data['RA_HC'][37], data['gauss_x_B3'][37]], [OW94_DEC, data['DEC_HC'][37], data['gauss_y_B3'][37]])
plt.grid()
plt.show()'''
