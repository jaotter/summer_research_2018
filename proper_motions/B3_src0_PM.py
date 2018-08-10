#SIMBAD name: HH508 [5,35,16.05,-5,23,7.2] and *tet01 Ori B [5,35,16.069,-5,23,6.96]

import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join, vstack
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits'))
obs = Table(ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt')) 

names = ['tet12', 'MAX', 'HC', 'B3', 'Fb']

ind = np.where(data['_idx_B3'] == 0)[0][0] 
ind_ref = np.where(data['_idx_B3'] == 3)[0][0]

RAs, DECs, dates, times_errs, ra_errs, dec_errs = get_info(ind, ind_ref, data, obs, names)
B3_ind = 3
ra_errs[B3_ind] = np.sqrt(data['x_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2)
dec_errs[B3_ind] = np.sqrt(data['y_err_B3'][ind]**2 + data['x_err_B3'][ind_ref]**2)

if data['PosFlag_MAX'][ind] == 0:
	MAX_poserr_ind = 0.1/3600
else:
	MAX_poserr_ind = 0.2/3600
if data['PosFlag_MAX'][ind_ref] == 0:
	MAX_poserr_ref = 0.1/3600
else:
	MAX_poserr_ref = 0.2/3600

MAX_poserr = np.sqrt(MAX_poserr_ind**2 + MAX_poserr_ref**2)
MAX_ind = 1
ra_errs[MAX_ind] = MAX_poserr
dec_errs[MAX_ind] = MAX_poserr

times = calc_times(dates)

v,pa = calc_pm(RAs, DECs, ra_errs, dec_errs, times, times_errs, 'HH 508', names)
print(v,pa)
plt.show()
