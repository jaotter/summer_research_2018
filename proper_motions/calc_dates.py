import numpy as np
def middle_date(st_mnth, st_yr, end_mnth, end_yr): #calculate the middle of the date range, and give the distance to each side (error)
	end_mnths = (end_mnth - st_mnth) + 12*(end_yr - st_yr)
	middle_mnths = end_mnths/2
	mnth = (middle_mnths+st_mnth) % 12
	yr = (middle_mnths+st_mnth)//12 + st_yr
	return (mnth, yr), middle_mnths/12

def calc_times(dates): #where dates is a list of tuples, (month, year)
	months = np.array([date[0] for date in dates])
	yrs = np.array([date[1] for date in dates])
	start = np.where(months+(yrs*12) == np.min(months+(yrs*12)))[0][0] #earliest date defined as t=0
	months = months - months[start]
	months[np.where(months < 0)] = months[np.where(months < 0)] + 12
	yrs = yrs - yrs[start]
	yrs[np.where(months < 0)] = yrs[np.where(months < 0)] - 1
	return yrs + months/12 #returns number of years after the earliest date
	

def sky_dist(x_i, y_i, x_f, y_f):
	del_x = x_f - x_i
	del_y = y_f - y_i
	dist = np.sqrt(del_x**2 + del_y**2)
	theta = (np.arctan(del_y/del_x)+np.pi/2)*(360/(2*np.pi))
	neg = np.where(del_x < 0)
	theta[neg] *= -1
	return dist, theta

def RA_to_deg(HH, MM, SS):
	return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
	return DD + (MM/60 + SS/(60*60))*np.sign(DD)

def get_info(ind, indref, data, obs_data, names):
	RAs = []
	DECs = []
	times = []
	times_errs = []
	ra_errs = []
	dec_errs = []
	for n in names:
		RAs.append(float(data['RA_'+n][ind] - data['RA_'+n][indref]))
		DECs.append(float(data['DEC_'+n][ind] - data['DEC_'+n][indref]))
		obs_ind = np.where(obs_data['name'] == n)
		times.append((float(obs_data['month'][obs_ind]), float(obs_data['year'][obs_ind])))
		times_errs.append(float(obs_data['t_err'][obs_ind].data[0]))
		ra_errs.append(np.sqrt(2)*float(obs_data['ra_err'][obs_ind].data[0])) #add factor of sqrt(2) bc adding in quadrature two sources with same error
		dec_errs.append(np.sqrt(2)*float(obs_data['dec_err'][obs_ind].data[0]))
	return RAs, DECs, times, times_errs, ra_errs, dec_errs

def get_errs(ind, ind_ref, data, obs_data, names, radio_only=True):
	times = []
	times_errs = []
	ra_errs = []
	dec_errs = []
	for n in names:
		obs_ind = np.where(obs_data['name'] == n)
		times.append((float(obs_data['month'][obs_ind]), float(obs_data['year'][obs_ind])))
		times_errs.append(float(obs_data['t_err'][obs_ind].data[0]))
		ra_errs.append(np.sqrt(2)*float(obs_data['ra_err'][obs_ind].data[0])) #add factor of sqrt(2) bc adding in quadrature two sources with same error
		dec_errs.append(np.sqrt(2)*float(obs_data['dec_err'][obs_ind].data[0]))
	if radio_only==True: #other data sets don't have the same naming convention for x/y errors or individual positional errors
		err_ind = np.where(np.isnan(ra_errs) == True)[0]
		for n_ind in err_ind:
			n = names[n_ind]
			ra_errs[n_ind] = np.sqrt(data['x_err_'+n][ind]**2 + np.sqrt(np.sum(data['x_err_'+n][ind_ref]**2))**2).quantity.value[0]
			dec_errs[n_ind] = np.sqrt(data['y_err_'+n][ind]**2 + np.sqrt(np.sum(data['y_err_'+n][ind_ref]**2))**2).quantity.value[0]
	return times, times_errs, ra_errs, dec_errs

