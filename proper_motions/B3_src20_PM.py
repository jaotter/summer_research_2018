#SIMBAD name [FW2011] J053514.656-052238.363, [5,35,14.656,-5,22,38.36]

import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits')) 
obs = Table(ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt')) 

ind = np.where(data['_idx_B3'] == 20)
ind_ref = np.where(data['_idx_B3'] == 3) #ref sources: 9(kinda bad), 97(good), 0(eh), 1(bad), 2(bad), 3(good)

names=['FW', 'HC', 'B3', 'LML', 'MLLA']

RAs, DECs, dates, times_errs, ra_errs, dec_errs = get_info(ind, ind_ref, data, obs, names)
nan_ind = np.where(np.isnan(ra_errs) == True)[0][0]
ra_errs[nan_ind] = np.sqrt(data['x_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2)
dec_errs[nan_ind] = np.sqrt(data['y_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2)

times = calc_times(dates)

v, pa = calc_pm(RAs, DECs, ra_errs, dec_errs, times, times_errs, 'src20', names)
print(v, pa)
plt.show()
