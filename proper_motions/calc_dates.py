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
	start = np.where(months+(yrs*12) == np.min(months+(yrs*12)))[0][0]
	months = months - months[start]
	months[np.where(months < 0)] = months[np.where(months < 0)] + 12
	yrs = yrs - yrs[start]
	yrs[np.where(months < 0)] = yrs[np.where(months < 0)] - 1
	return yrs + months/12
	

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
