import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from functools import reduce
from calc_dates import *
from os import sys

def f(B, x): #m is list of parameters (slope and intercept)
	return B[0]*x + B[1]

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits'))
obs = Table(ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt')) 

names = np.array(['B3', 'B6', 'Fb']) #'B3', 'B6', 'A340', 'A470', 'Fb'

ind_list = []
for n in names:
	ind_list.append(np.where(np.isnan(data['RA_'+n]) == False)[0])
all_ind = reduce(np.intersect1d, ind_list) 
all_ind = np.delete(all_ind, np.where(all_ind == 5)[0])

rem_list = [np.where(all_ind==10)[0][0], np.where(all_ind==20)[0][0], np.where(all_ind==54)[0][0]] #np.where(all_ind==5)[0][0], #src n, src I, BN, all of which are known to have high PM
#rem_list = [np.where(all_ind==10), np.where(all_ind==20)]
ref_ind = np.delete(all_ind, rem_list) #reference sources - detected in all catalogs and no high proper motions
arr = np.arange(len(all_ind))
ref_all_ind = np.delete(arr, rem_list) #reference sources in all_ind

dates, times_errs, ra_errs, dec_errs = get_errs(all_ind, ref_ind, data, obs, names)
times = calc_times(dates)

coord_table = Table()
coord_table['D_ID'] = data['D_ID'][all_ind]

for i,n in enumerate(names): #make calibration position
	#first average ra/dec of reference sources
	ra_ref = np.mean(data['RA_'+n][ref_ind])
	de_ref = np.mean(data['DEC_'+n][ref_ind])
	#now subtract avg from remaining sources
	RAs = data['RA_'+n][all_ind] - ra_ref
	DECs = data['DEC_'+n][all_ind] - de_ref
	coord_table['RA_'+n] = RAs
	coord_table['DEC_'+n] = DECs
	try:
		coord_table['RA_err_'+n] = np.sqrt(data['x_err_'+n][all_ind]**2 + np.sum(data['x_err_'+n][ref_ind]**2))
		coord_table['DEC_err_'+n] = np.sqrt(data['y_err_'+n][all_ind]**2 + np.sum(data['y_err_'+n][ref_ind]**2))
	except KeyError:
		coord_table['RA_err_'+n] = np.full(len(coord_table), ra_errs[i]*np.sqrt(len(ref_ind)+1))
		coord_table['DEC_err_'+n] = np.full(len(coord_table), dec_errs[i]*np.sqrt(len(ref_ind)+1))

#B3/B6 at same time, so avg positions
coord_table['RA_B36'] = (coord_table['RA_B3'] + coord_table['RA_B6'])/2
coord_table['DEC_B36'] = (coord_table['DEC_B3'] + coord_table['DEC_B6'])/2
coord_table['RA_err_B36'] = np.sqrt(coord_table['RA_err_B3']**2 + coord_table['RA_err_B6']**2)
coord_table['DEC_err_B36'] = np.sqrt(coord_table['DEC_err_B3']**2 + coord_table['DEC_err_B6']**2)

time = 5*u.yr #5 years 0 months

ra_pm = ((coord_table['RA_B36'] - coord_table['RA_Fb'])*u.degree).to(u.mas)/time
ra_pm = ra_pm*np.cos(0.25*data['DEC_B3'][all_ind]+0.25*data['DEC_B6'][all_ind]+0.5*data['DEC_Fb'][all_ind])
dec_pm = ((coord_table['DEC_B36'] - coord_table['DEC_Fb'])*u.degree).to(u.mas)/time
ra_pm_err = np.std(ra_pm[ref_all_ind]) #only use reference source in std dev
dec_pm_err = np.std(dec_pm[ref_all_ind])

tot_pm = np.sqrt(ra_pm**2 + dec_pm**2)
tot_pm_err = np.sqrt(ra_pm_err**2 + dec_pm_err**2)

pa = np.pi/2*u.rad - np.arctan((dec_pm/ra_pm).decompose())
pa[np.where(ra_pm < 0)] = -np.pi*u.rad+pa[np.where(ra_pm < 0)]
#tab = Table((coord_table['D_ID'], tot_pm, tot_pm_err), names=('D_ID', 'pm', 'pm_err'))

plt.errorbar(coord_table['RA_B36'][rem_list], coord_table['DEC_B36'][rem_list], xerr=coord_table['RA_err_B36'][rem_list], yerr= coord_table['DEC_err_B36'][rem_list], linestyle='',label='t=5')
plt.errorbar(coord_table['RA_Fb'][rem_list], coord_table['DEC_Fb'][rem_list], xerr=coord_table['RA_err_Fb'][rem_list], yerr=coord_table['DEC_err_Fb'][rem_list], linestyle='',label='t=0')

plt.text(coord_table['RA_B36'][rem_list[0]], coord_table['DEC_B36'][rem_list[0]], 'n')
plt.text(coord_table['RA_B36'][rem_list[1]], coord_table['DEC_B36'][rem_list[1]], 'I')
plt.text(coord_table['RA_B36'][rem_list[2]], coord_table['DEC_B36'][rem_list[2]], 'BN')

plt.xlabel('offset in RA (deg)')
plt.ylabel('offset in DEC (deg)')

plt.xlim(0, -.002)

plt.legend()


#Quiver plot of proper motions
fl = fits.open('/lustre/aoc/students/jotter/directory/Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits')
header = fl[0].header
img = fl[0].data
mywcs = WCS(header).celestial
#img = img[::-1]
ax = plt.subplot(projection=mywcs)
ax.imshow(img, vmin=0.5e-4, vmax=1e-3, origin='lower')

#add measurements using multi wavelength data
alt_ind = np.array([np.where(data['_idx_B3'] == 20)[0][0], np.where(data['_idx_B3'] == 62)[0][0], np.where(data['_idx_B3'] == 105)[0][0]])
new_ind = np.concatenate((all_ind, alt_ind))


alt_pm = [8.11, 5.40, 7.57]
alt_pa = np.array([-1.28, 2.63, 1.81])*u.rad
alt_pmx = np.sin(alt_pa)*alt_pm*np.cos(data['DEC_B3'][alt_ind])
alt_pmy = np.cos(alt_pa)*alt_pm

pa = np.concatenate((pa, alt_pa))
pa = pa.value*u.rad
tot_pm = np.concatenate((tot_pm, alt_pm))

ra_pm = np.concatenate((ra_pm, alt_pmx))
dec_pm = np.concatenate((dec_pm, alt_pmy))

arrow_x = data['RA_B3'][new_ind]
arrow_y = data['DEC_B3'][new_ind]

arrow_x_pix, arrow_y_pix = mywcs.all_world2pix(arrow_x, arrow_y, 1)

arrow_xlen = -1*np.sin(pa)*tot_pm*np.cos(data['DEC_B3'][new_ind])
arrow_ylen = np.cos(pa)*tot_pm

src_names = np.array([str(i) for i in data['D_ID'][new_ind].data])
src_names[rem_list] = ['I', 'BN', '54'] #'n'

for i in range(len(src_names)):
	plt.text(arrow_x_pix[i]-5, arrow_y_pix[i]+10, src_names[i], color='r')
Q = plt.quiver(arrow_x_pix, arrow_y_pix, arrow_xlen, arrow_ylen, color='g')
plt.quiverkey(Q, 0.1, 0.1, 10, '10 mas/yr', labelcolor='r')

ax.set_xlabel('RA')
ax.set_ylabel('DEC')
plt.show()

data['mu_alpha_cos_dec'] = np.full(len(data), np.nan)
data['mu_alpha_err'] = np.full(len(data), np.nan)
data['mu_delta'] = np.full(len(data), np.nan)
data['mu_delta_err'] = np.full(len(data), np.nan)
data['pm_pa'] = np.full(len(data), np.nan)

data['mu_alpha_cos_dec'][new_ind] = arrow_xlen
data['mu_alpha_err'][new_ind] = ra_pm_err #NOTE: INCORRECT ERRORS FOR SOURCES MEASURED WITH MULTI-WAVELENGTH CATALOGS POSITION FITTING
data['mu_delta'][new_ind] = arrow_ylen
data['mu_delta_err'][new_ind] = dec_pm_err
data['pm_pa'][new_ind] = pa

data.write('/lustre/aoc/students/jotter/dendro_catalogs/mm_IR_data_pm.fits', overwrite=True)

