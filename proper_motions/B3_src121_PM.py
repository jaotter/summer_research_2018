#SIMBAD name: [MRD2012] 2308

import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits')) 
obs = Table(ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt')) 

names = ['OW94', 'HC', 'B3', 'MRD', 'COUP']

ind = np.where(data['_idx_B3'] == 121)
ind_ref = np.where(data['_idx_B3'] == 82)

RAs, DECs, dates, times_errs, ra_errs, dec_errs = get_info(ind, ind_ref, data, obs, names)
COUP_ind = 4
ra_errs[COUP_ind] = np.sqrt(data['PosErr_COUP'][ind]**2 + data['PosErr_COUP'][ind_ref]**2).quantity.value[0]/3600
dec_errs[COUP_ind] = np.sqrt(data['PosErr_COUP'][ind]**2 + data['PosErr_COUP'][ind_ref]**2).quantity.value[0]/3600
B3_ind = 2
ra_errs[B3_ind] = np.sqrt(data['x_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2).quantity.value
dec_errs[B3_ind] = np.sqrt(data['y_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2).quantity.value

times = calc_times(dates)

v,pa = calc_pm(RAs, DECs, ra_errs, dec_errs, times, times_errs, 'src121')
print(v,pa)
plt.show()
