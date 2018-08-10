#SIMBAD name: [SCE2006] 15 [5,35,14.656,-5,22,38.36]

import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits')) 
obs = Table(ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt')) 

names=['COUP', 'HC', 'B3', 'FDM', 'MLLA']

ind = np.where(data['_idx_B3'] == 62)
ind_ref = np.where(data['_idx_B3'] == 3) #ref sources: 32(meh),82(bad),94(good),0(bad),1(meh),3(ok)

RAs, DECs, dates, times_errs, ra_errs, dec_errs = get_info(ind, ind_ref, data, obs, names)

B3_ind = 2
ra_errs[B3_ind] = np.sqrt(data['x_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2).quantity.value[0]
dec_errs[B3_ind] = np.sqrt(data['y_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2).quantity.value[0]
COUP_ind = 0
ra_errs[COUP_ind] = np.sqrt(data['PosErr_COUP'][ind]**2 + data['PosErr_COUP'][ind_ref]**2).quantity.value[0]/3600
dec_errs[COUP_ind] = np.sqrt(data['PosErr_COUP'][ind]**2 + data['PosErr_COUP'][ind_ref]**2).quantity.value[0]/3600
FDM_ind = 3
ra_errs[FDM_ind] = np.sqrt(data['PosErr_FDM'][ind]**2 + data['PosErr_FDM'][ind_ref]**2).quantity.value[0]/3600
dec_errs[FDM_ind] = np.sqrt(data['PosErr_FDM'][ind]**2 + data['PosErr_FDM'][ind_ref]**2).quantity.value[0]/3600

times = calc_times(dates)

v,pa = calc_pm(RAs, DECs, ra_errs, dec_errs, times, times_errs, 'src62', names)
print(v,pa)
plt.show()

