import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import numpy as np


#data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/allbands_catalog.fits')
data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'))
fig = plt.figure()
for row in range(len(data)):
	x_470GHz = data['gauss_x_470GHz'][row]
	y_470GHz = data['gauss_y_470GHz'][row]
	x_340GHz = data['gauss_x_340GHz'][row]
	y_340GHz = data['gauss_y_340GHz'][row]
	x_B3 = data['gauss_x_B3'][row] 
	y_B3 = data['gauss_y_B3'][row] 
	x_B6 = data['gauss_x_B6'][row]
	y_B6 = data['gauss_y_B6'][row]  
	RA_HC = data['RA_HC'][row]
	DEC_HC = data['DEC_HC'][row]
	RA_RRS = data['RA_RRS'][row]
	DEC_RRS = data['DEC_RRS'][row]
	plot_ra = np.array([RA_HC, RA_RRS, x_470GHz, x_340GHz, x_B3, x_B6])
	plot_ra = plot_ra[np.where(np.isnan(plot_ra) == False)]
	plot_dec = np.array([DEC_HC, DEC_RRS, y_470GHz, y_340GHz, y_B3, y_B6])
	plot_dec = plot_dec[np.where(np.isnan(plot_dec) == False)]
	plt.plot(plot_ra, plot_dec, color='r')


plt.scatter(data['gauss_x_340GHz'], data['gauss_y_340GHz'], c='b', alpha=0.4, label='340GHz')
plt.scatter(data['gauss_x_470GHz'], data['gauss_y_470GHz'], c='g', alpha=0.4, label='470GHz')
plt.scatter(data['gauss_x_B3'], data['gauss_y_B3'], c='r', alpha=0.4, label='B3')
plt.scatter(data['gauss_x_B6'], data['gauss_y_B6'], c='y', alpha=0.4, label='B6')
plt.scatter(data['RA_HC'], data['DEC_HC'], c='violet', alpha=0.4, label='HC IR')
plt.scatter(data['RA_RRS'], data['DEC_RRS'], c='k', alpha=0.4, label='RRS IR')

plt.xlabel('RA (deg)')
plt.ylabel('DEC (deg)')
plt.legend()
plt.show()
