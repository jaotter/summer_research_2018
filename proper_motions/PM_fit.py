from scipy.odr import *
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from calc_dates import *

def f(B, x): #m is list of parameters (slope and intercept)
	return B[0]*x + B[1]

def calc_pm(RAs, DECs, RA_err, DEC_err, times, times_err, src_name, names, slope_guess=1, int_guess=-100, published_pm = None):
	if len(RAs) == 2: #Currently having some problems with scipy.odr with only two elements
		time = np.abs(times[1] - times[0])
		ra_pm = (RAs[0] - RAs[1])/time
		ra_pm = ra_pm*np.cos(0.5*DECs[0]+0.5*DECs[1])
		dec_pm = (DECs[0] - DECs[1])/time

		tot_pm = np.sqrt(ra_pm**2 + dec_pm**2)*u.deg.to(u.mas)

		pa = (np.pi/2) - np.arctan(dec_pm/ra_pm)
		if ra_pm < 0:
			pa = -np.pi + pa
		return tot_pm, pa


	src_name = str(src_name)
	linear = Model(f)
	ind = np.where(np.abs(RAs) < 100000)[0] #in the table, masked RA/DEC values are 9999999, so filter those out
	RAs = np.array(RAs)[ind]
	DECs = np.array(DECs)[ind]
	RA_err = np.array(RA_err)[ind]
	DEC_err = np.array(DEC_err)[ind]
	times = np.array(times)[ind]
	times_err = np.array(times_err)[ind]
	names = np.array(names)[ind]
	mydata = RealData(RAs, DECs, sx=RA_err, sy=DEC_err)
	myodr = ODR(mydata, linear, beta0=[1,-100])
	myoutput=myodr.run()

	slope = myoutput.beta[0]
	inter = myoutput.beta[1]
	
	time_zero_ind = np.where(times == np.min(times))
	time_f_ind = np.where(times == np.max(times))
	pa = np.pi/2 - np.arctan(slope)
	if RAs[time_f_ind] < RAs[time_zero_ind]:
		pa = -np.pi + pa

	dists = [] #distance of point from fitted line
	line_xy = [] #point along fitted line

	for i in range(len(RAs)):
		dist = np.abs(-slope*RAs[i] + DECs[i] - inter)/np.sqrt(slope**2 + 1)
		dists.append(dist)
		
		x_line = (RAs[i] + slope*DECs[i] - slope*inter)/(1 + slope**2)
		y_line = (slope*RAs[i] + DECs[i]*(slope**2) + inter)/(1 + slope**2)
		line_xy.append((x_line, y_line))
	
	dists=(dists*u.degree).to(u.mas)

	min_ind = np.where(times == np.min(times))[0][0]
	
	pm_dist = []
	for ind, pt in enumerate(line_xy):
		d = np.sqrt((line_xy[min_ind][0] - pt[0])**2 + (line_xy[min_ind][1] - pt[1])**2) #distance along vector i.e. proper motion distance
		t = times[ind] - times[min_ind]
		pm_dist.append(d)

	vel_data = RealData(times, pm_dist.value, sx=times_err, sy=dists.value)	
	vel_odr = ODR(vel_data, linear, beta0=[1,1])
	vel_out = vel_odr.run()
	
	pm = vel_out.beta[0]*u.deg.to(u.mas)
	pm_err = vel_out.sd_beta[0]
	x_vals2 = np.linspace(np.min(times), np.max(times), 10)

	

	plt.figure(figsize=(8,6))	
	plt.errorbar(RAs, DECs, xerr=RA_err, yerr=DEC_err, linestyle='', marker='o')
	plt.xlabel(r'$\Delta$ RA (degree)')
	plt.ylabel(r'$\Delta$ DEC (degree)')
	x_vals = np.linspace(np.nanmin(RAs), np.nanmax(RAs),10)
	plt.plot(x_vals, f(myoutput.beta, x_vals))
	for i in range(len(RAs)):
		plt.text(RAs[i], DECs[i], names[i])
	if published_pm is not None:
		mag = published_pm[0]
		p_pa = (published_pm[1]*u.degree).to(u.rad)
		plt.plot(x_vals, (1/np.tan(p_pa))*x_vals + (myoutput.beta[0]*x_vals[5] + myoutput.beta[1] - (1/np.tan(p_pa))*x_vals[5]), label='published PM, pa = {:.1f}'.format(published_pm[1]))
	plt.legend()
	plt.savefig('/users/jotter/summer_research_2018/proper_motions/PM_plots/'+src_name+'_ra-dec.png',dpi=300)

	plt.figure(figsize=(8,6))
	plt.plot(x_vals2, f(vel_out.beta, x_vals2), label=r'$\mu = {:.1f} \pm {:f}$'.format(pm, pm_err))
	plt.errorbar(times, pm_dist.value, xerr=times_err, yerr=dists.value, linestyle='', marker='o')	
	plt.xlabel('time (years)')
	plt.ylabel('source movement (mas)')
	plt.xlim(-0.5, np.nanmax(times)+1)
	plt.ylim(-0.1*np.nanmax(pm_dist.value), 1.1*np.nanmax(pm_dist.value))
	for i in range(len(times)):
		plt.text(times[i], pm_dist.value[i], names[i])
	if published_pm is not None:
		plt.plot(x_vals2, mag*x_vals2, label='published PM, '+r'$\mu = {:.1f}$'.format(mag))
	plt.legend()

	plt.savefig('/users/jotter/summer_research_2018/proper_motions/PM_plots/'+src_name+'_time_dist.png',dpi=300)

	return vel_out.beta[0], pa






