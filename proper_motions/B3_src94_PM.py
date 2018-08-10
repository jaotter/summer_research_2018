#SIMBAD name: [MRD2012] 2297

import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits')) 
obs = Table(ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt')) 

names=['COUP', 'HC', 'B3', 'MLLA', 'OW94', 'MAX', 'MRD', 'RRS', 'Fb']

ind = np.where(data['_idx_B3'] == 94)
ind_ref = np.where(data['_idx_B3'] == 3) #ref sources:0(bad),8(bad),7(kinda bad),82(okay),97(ok),3

RAs, DECs, dates, times_errs, ra_errs, dec_errs = get_info(ind, ind_ref, data, obs, names)
COUP_ind = 0
ra_errs[COUP_ind] = np.sqrt(data['PosErr_COUP'][ind]**2 + data['PosErr_COUP'][ind_ref]**2).quantity.value[0]/3600
dec_errs[COUP_ind] = np.sqrt(data['PosErr_COUP'][ind]**2 + data['PosErr_COUP'][ind_ref]**2).quantity.value[0]/3600
B3_ind = 2
ra_errs[B3_ind] = np.sqrt(data['x_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2).quantity.value
dec_errs[B3_ind] = np.sqrt(data['y_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2).quantity.value
if data['PosFlag_MAX'][ind] == 0:
	MAX_poserr_ind = 0.1/3600
else:
	MAX_poserr_ind = 0.2/3600
if data['PosFlag_MAX'][ind_ref] == 0:
	MAX_poserr_ref = 0.1/3600
else:
	MAX_poserr_ref = 0.2/3600
MAX_poserr = np.sqrt(MAX_poserr_ind**2 + MAX_poserr_ref**2)
MAX_ind = 5
ra_errs[MAX_ind] = MAX_poserr
dec_errs[MAX_ind] = MAX_poserr
times = calc_times(dates)

v,pa = calc_pm(RAs, DECs, ra_errs, dec_errs, times, times_errs, 'src94', names)
print(v,pa)
plt.show()

