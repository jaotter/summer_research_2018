import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import numpy as np


#data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/allbands_catalog.fits')
data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'))
IR_data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_only_catalog.fits'))
fig = plt.figure()
#ref_ind = np.where(data['_idx_B3']==89) #index of BN, one of only sources in IR and radio catalogs
ref_ind = np.where(data['_idx_B3']==40) 
x_470GHz = data['gauss_x_470GHz'] - data['gauss_x_470GHz'][ref_ind]
y_470GHz = data['gauss_y_470GHz'] - data['gauss_y_470GHz'][ref_ind]
x_340GHz = data['gauss_x_340GHz'] - data['gauss_x_340GHz'][ref_ind]
y_340GHz = data['gauss_y_340GHz'] - data['gauss_y_340GHz'][ref_ind]
x_B3 = data['gauss_x_B3'] - data['gauss_x_B3'][ref_ind]
y_B3 = data['gauss_y_B3'] - data['gauss_y_B3'][ref_ind] 
x_B6 = data['gauss_x_B6'] - data['gauss_x_B6'][ref_ind]
y_B6 = data['gauss_y_B6'] - data['gauss_y_B6'][ref_ind]
RA_HC = data['RA_HC'] - data['RA_HC'][ref_ind]
DEC_HC = data['DEC_HC'] - data['DEC_HC'][ref_ind]
RA_RRS = data['RA_RRS'] - data['RA_RRS'][ref_ind]
DEC_RRS = data['DEC_RRS'] - data['DEC_RRS'][ref_ind]
for row in range(len(data)):
	#plot_ra = np.array([RA_HC[row], RA_RRS[row], x_470GHz[row], x_340GHz[row], x_B6[row], x_B3[row]])
	plot_ra = np.array([RA_HC[row], RA_RRS[row], x_470GHz[row], x_340GHz[row], x_B6[row]])
	plot_ra = plot_ra[np.where(np.isnan(plot_ra) == False)]
	x_err = np.array([data['RA_HC'][row] - data['RA_HC'][row], data['RA_RRS'][row] - data['RA_RRS'][row], data['x_err_470GHz'][row], data['x_err_340GHz'][row], data['x_err_B6'][row]])
	x_err = x_err[np.where(np.isnan(plot_ra) == False)]
	#plot_dec = np.array([DEC_HC[row], DEC_RRS[row], y_470GHz[row], y_340GHz[row], y_B6[row], y_B3[row]])
	plot_dec = np.array([DEC_HC[row], DEC_RRS[row], y_470GHz[row], y_340GHz[row], y_B6[row]])
	plot_dec = plot_dec[np.where(np.isnan(plot_dec) == False)]
	y_err = np.array([data['RA_HC'][row] - data['RA_HC'][row], data['RA_RRS'][row] - data['RA_RRS'][row], data['y_err_470GHz'][row], data['y_err_340GHz'][row], data['y_err_B6'][row]])
	y_err = y_err[np.where(np.isnan(plot_dec) == False)]
	#plt.plot(plot_ra, plot_dec, color='r')
	#plt.plot(plot_ra - x_B3[row], plot_dec - y_B3[row], marker='*', linestyle='-')
	#plt.errorbar(plot_ra - x_B3[row], plot_dec - y_B3[row], xerr=x_err, yerr=y_err, marker='*')


plt.errorbar(data['gauss_x_340GHz'] - data['gauss_x_340GHz'][ref_ind] - x_B3, data['gauss_y_340GHz'] - data['gauss_y_340GHz'][ref_ind] - y_B3, xerr=data['x_err_340GHz'], yerr=data['y_err_340GHz'], c='b', alpha=0.4, label='340GHz', linestyle='', marker='*')
plt.errorbar(data['gauss_x_470GHz'] - data['gauss_x_470GHz'][ref_ind] - x_B3, data['gauss_y_470GHz'] - data['gauss_y_470GHz'][ref_ind] - y_B3, xerr=data['x_err_470GHz'], yerr=data['y_err_470GHz'], c='g', alpha=0.4, label='470GHz', linestyle='', marker='*')
#plt.errorbar(data['gauss_x_B3'] - data['gauss_x_B3'][ref_ind], data['gauss_y_B3'] - data['gauss_y_B3'][ref_ind], xerr=data['x_err_B3'], yerr=data['y_err_B3'], c='r', alpha=0.4, label='B3', linestyle='', marker='*')
plt.errorbar(data['gauss_x_B6'] - data['gauss_x_B6'][ref_ind] - x_B3, data['gauss_y_B6'] - data['gauss_y_B6'][ref_ind] - y_B3, xerr=data['x_err_B6'], yerr=data['y_err_B6'], c='y', alpha=0.4, label='B6', linestyle='', marker='*')
plt.scatter(data['RA_HC'] - data['RA_HC'][ref_ind] - x_B3, data['DEC_HC'] - data['DEC_HC'][ref_ind] - y_B3, c='violet', alpha=0.4, label='HC IR')
plt.scatter(data['RA_RRS'] - data['RA_RRS'][ref_ind] - x_B3, data['DEC_RRS'] - data['DEC_RRS'][ref_ind] - y_B3, c='k', alpha=0.4, label='RRS IR')

#plt.scatter(IR_data['RA_HC'], IR_data['DEC_HC'], c='violet', alpha=0.4, label='HC IR')
#plt.scatter(IR_data['RA_RRS'], IR_data['DEC_RRS'], c='k', alpha=0.4, label='RRS IR')

plt.xlabel('Delta RA (deg)')
plt.ylabel('Delta DEC (deg)')
plt.legend()
plt.grid()
plt.show()
