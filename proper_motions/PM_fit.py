from scipy.odr import *
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from calc_dates import *

def f(B, x): #m is list of parameters (slope and intercept)
	return B[0]*x + B[1]

def calc_pm(RAs, DECs, RA_err, DEC_err, times, times_err):
	linear = Model(f)
	ind = np.where(np.isnan(RAs) == False)[0]
	RAs = np.array(RAs)[ind]
	DECs = np.array(DECs)[ind]
	RA_err = np.array(RA_err)[ind]
	DEC_err = np.array(DEC_err)[ind]
	times = np.array(times)[ind]
	times_err = np.array(times_err)[ind]
	mydata = RealData(RAs, DECs, sx=RA_err, sy=DEC_err)
	myodr = ODR(mydata, linear, beta0=[1,-100])
	myoutput=myodr.run()
	myoutput.pprint()

	slope = myoutput.beta[0]
	inter = myoutput.beta[1]
	
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
	pm_dist = (pm_dist*u.degree).to(u.mas)
	vel_data = RealData(times, pm_dist.value, sx=times_err, sy=dists.value)	
	vel_odr = ODR(vel_data, linear, beta0=[1,1])
	vel_out = vel_odr.run()
	
	pm = vel_out.beta[0]
	pm_err = vel_out.sd_beta[0]
	x_vals2 = np.linspace(np.min(times), np.max(times), 10)

	

	plt.figure()	
	plt.errorbar(RAs, DECs, xerr=RA_err, yerr=DEC_err, linestyle='', marker='o')
	plt.xlabel(r'$\Delta$ RA (degree)')
	plt.ylabel(r'$\Delta$ DEC (degree)')
	x_vals = np.linspace(np.nanmin(RAs), np.nanmax(RAs),10)
	plt.plot(x_vals, f(myoutput.beta, x_vals))

	plt.figure()
	plt.plot(x_vals2, f(vel_out.beta, x_vals2), label=r'$\mu = {:.1f} \pm {:f}$'.format(pm, pm_err))
	plt.errorbar(times, pm_dist.value, xerr=times_err, yerr=dists.value, linestyle='', marker='o')	
	plt.xlabel('time (years)')
	plt.ylabel('source movement (degree)')
	plt.xlim(-0.5, np.nanmax(times)+1)
	plt.ylim(-0.1*np.nanmax(pm_dist.value), 1.1*np.nanmax(pm_dist.value))
	plt.legend()
	
	return vel_out.beta[0]






