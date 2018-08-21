import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from functools import reduce
from calc_dates import *
from PM_fit import *
from os import sys

#data = '/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits'
#obs_data = '/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt'
#names = np.array(['B3', 'B6', 'HC', 'Fb']) #'B3', 'B6', 'A340', 'A470', 'Fb'
#rem_list = [5,10,20]

def calc_pms(data, obs_data, names, rem_list):
	#parameters:
	#data - master catalog of all matched catalogs, with RA/DEC of each source
	#obs_data - master catalog of time of each observation, error in time, and error in ra/dec
	#names - names of catalogs to use
	#rem_list - sources with known high PM to exclude from the reference position (orion: 5 - n, 10 - I, 20 - BN)



	data = Table.read(data)
	obs = Table(ascii.read(obs_data))

	ind_list = [] #all_ind is the list of indices of sources common to all catalogs in names
	for n in names:
		ind_list.append(np.where(np.isnan(data['RA_'+n]) == False)[0])
	all_ind = reduce(np.intersect1d, ind_list) 

	rm_inds = []
	for rl in rem_list:
		rlind = np.where(all_ind==rl)[0]
		if len(rlind) > 0:
			rm_inds.append(rlind)
	ref_ind = np.delete(all_ind, rm_inds) #reference sources - detected in all catalogs and no high proper motions
	arr = np.arange(len(all_ind))
	ref_all_ind = np.delete(arr, rm_inds) #reference source indices in all_ind

	dates, times_errs, ra_errs, dec_errs = get_errs(all_ind, ref_ind, data, obs, names) #function from calc_dates.py which gives the observation info (date, date error, ra/dec error)
	times = calc_times(dates) #function which turns dates into times, where the earliest date is t=0. Units are months

	coord_table = Table() 
	coord_table['D_ID'] = data['D_ID'][all_ind]

	for i,n in enumerate(names): #make calibration position by averaging RA/DEC of reference sources - then subtract from the other sources
		#first average ra/dec of reference sources
		ra_ref = np.mean(data['RA_'+n][ref_ind])
		de_ref = np.mean(data['DEC_'+n][ref_ind])
		#now subtract avg from remaining sources
		RAs = data['RA_'+n][all_ind] - ra_ref
		DECs = data['DEC_'+n][all_ind] - de_ref
		coord_table['RA_'+n] = RAs
		coord_table['DEC_'+n] = DECs
		try: #this checks if there are source specific position errors in data, otherwise use the errors from obs_data
			coord_table['RA_err_'+n] = np.sqrt(data['x_err_'+n][all_ind]**2 + np.sum(data['x_err_'+n][ref_ind]**2))
			coord_table['DEC_err_'+n] = np.sqrt(data['y_err_'+n][all_ind]**2 + np.sum(data['y_err_'+n][ref_ind]**2))
		except KeyError:
			coord_table['RA_err_'+n] = np.full(len(coord_table), ra_errs[i]*np.sqrt(len(ref_ind)+1))
			coord_table['DEC_err_'+n] = np.full(len(coord_table), dec_errs[i]*np.sqrt(len(ref_ind)+1))

	#measurements at same time should have averaged positions

	same_epoch = [] #list of list of indices where measurements were taken at the same time
	un, counts = np.unique(times, return_counts=True)
	for n in range(len(counts)):
		if counts[n] > 1:
			same_epoch.append(np.where(times == un[n]))

	full_names = names #full_names will include the concatenated names for same-time observations
	full_times = times
	full_times_err = times_errs
	for se in same_epoch[0]: #now go through same time measurements and average positions
		nm = ''
		ras = []
		decs = []
		ra_errs = []
		dec_errs = []
		for i in se:
			nm = nm+str(names[i])
			ras.append(coord_table['RA_'+names[i]])
			decs.append(coord_table['DEC_'+names[i]])
			ra_errs.append(coord_table['RA_err_'+names[i]])
			dec_errs.append(coord_table['DEC_err_'+names[i]])
			names_ind = np.where(full_names == names[i])[0]
			full_names = np.delete(full_names, names_ind) #remove specific element not index
			full_times = np.delete(full_times, names_ind)
			full_times_err = np.delete(full_times_err, names_ind)
		coord_table['RA_'+nm] = [np.mean(np.array(ras)[:,row]) for row in range(len(ras[0]))] #average of positions
		coord_table['DEC_'+nm] = [np.mean(np.array(decs)[:,row]) for row in range(len(decs[0]))]
		coord_table['RA_err_'+nm] = [np.sqrt(np.sum(np.array(ra_errs)[:,row]**2)) for row in range(len(ra_errs[0]))] #add in quadrature
		coord_table['DEC_err_'+nm] = [np.sqrt(np.sum(np.array(dec_errs)[:,row]**2)) for row in range(len(dec_errs[0]))]
		full_names = np.append(full_names, nm)
		full_times = np.append(full_times, times[se[0]])
		full_times_err = np.append(full_times_err, times_errs[se[0]])

	tot_pm_arr = []
	pa_arr = []
	for src in range(len(coord_table)):
		RAs = []
		DECs = []
		RA_err = []
		DEC_err = []

		for name in full_names:
			RAs.append(coord_table['RA_'+name][src])
			DECs.append(coord_table['DEC_'+name][src])
			RA_err.append(coord_table['RA_err_'+name][src])
			DEC_err.append(coord_table['DEC_err_'+name][src])
			names_ind = np.where(names == name)[0]
		
		pm, pa = calc_pm(RAs, DECs, RA_err, DEC_err, full_times, full_times_err, coord_table['D_ID'][src], full_names)
		tot_pm_arr.append(pm)
		pa = pa*u.rad.to(u.deg)
		pa_arr.append(pa)

	coord_table['pm'] = tot_pm_arr
	coord_table['pa'] = pa_arr



