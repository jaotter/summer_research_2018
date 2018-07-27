import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle

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

time = {'HC':0, 'RRS':8, '340GHz':17, '470GHz':15, 'B3':18, 'B6':18, 'B7':18, 'Fb':6}

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'))
ref_ind = np.where(data['_idx_B3']==32) 

x_470GHz = data['gauss_x_470GHz'] - data['gauss_x_470GHz'][ref_ind]
y_470GHz = data['gauss_y_470GHz'] - data['gauss_y_470GHz'][ref_ind]
x_340GHz = data['gauss_x_340GHz'] - data['gauss_x_340GHz'][ref_ind]
y_340GHz = data['gauss_y_340GHz'] - data['gauss_y_340GHz'][ref_ind]
x_B3 = data['gauss_x_B3'] - data['gauss_x_B3'][ref_ind]
y_B3 = data['gauss_y_B3'] - data['gauss_y_B3'][ref_ind] 
x_B6 = data['gauss_x_B6'] - data['gauss_x_B6'][ref_ind]
y_B6 = data['gauss_y_B6'] - data['gauss_y_B6'][ref_ind]
RA_Fb = data['RA_Fb'] - data['RA_Fb'][ref_ind]
DEC_Fb = data['DEC_Fb'] - data['DEC_Fb'][ref_ind]
RA_HC = data['RA_HC'] - data['RA_HC'][ref_ind]
DEC_HC = data['DEC_HC'] - data['DEC_HC'][ref_ind]
RA_RRS = data['RA_RRS'] - data['RA_RRS'][ref_ind]
DEC_RRS = data['DEC_RRS'] - data['DEC_RRS'][ref_ind]

dist_table = Table()#names=('D_ID', 'v_HC', 'pa_HC', 'v_RRS', 'pa_RRS', 'v_340GHz', 'pa_340GHz', 'v_470GHz', 'pa_470GHz', 'd_B6', 'pa_B6', 'd_B7', 'pa_B7'))

#for row in range(len(data)):
x_vals = [RA_RRS, RA_HC, x_470GHz, x_340GHz, x_B6, RA_Fb]
y_vals = [DEC_RRS, DEC_HC, y_470GHz, y_340GHz, y_B6, DEC_Fb]
names = ['RRS', 'HC', '470GHz', '340GHz', 'B6', 'Fb']
HC_ra_err = (0.017*u.arcsecond).to(u.degree)
HC_dec_err = (0.028*u.arcsecond).to(u.degree)
x_err = [0, np.full(len(data), HC_ra_err), data['x_err_470GHz'], data['x_err_340GHz'], data['x_err_B6'], data['e_RAs_Fb']]
y_err = [0, np.full(len(data), HC_dec_err), data['y_err_470GHz'], data['y_err_340GHz'], data['y_err_B6'], data['e_DEs_Fb']]
for i in range(len(x_vals)):
	#ind = np.where(np.isnan(x_vals[i]) == False)
	dist_arr, theta_arr = sky_dist(x_vals[i], y_vals[i], x_B3, y_B3)
	dist_err_ra = np.sqrt(x_err[i]**2 + data['x_err_B3']**2)
	dist_err_dec = np.sqrt(y_err[i]**2 + data['y_err_B3']**2)
	delta_t = time['B3'] - time[names[i]]
	if delta_t != 0:
		v = dist_arr/delta_t * u.degree / u.year
		v = v.to(u.mas/u.year)
		v_err_x = dist_err_ra/delta_t * u.degree/u.year
		v_err_x = v_err_x.to(u.mas/u.year)
		v_err_y = dist_err_dec/delta_t * u.degree/u.year
		v_err_y = v_err_y.to(u.mas/u.year)
		dist_table['v_'+names[i]] = v
		dist_table['v_err_ra_'+names[i]] = v_err_x
		dist_table['v_err_dec_'+names[i]] = v_err_y
	else:
		dist_table['d_'+names[i]] = dist_arr
	dist_table['pa_'+names[i]] = theta_arr
dist_table['D_ID'] = data['D_ID']



img, header = fits.getdata('/lustre/aoc/students/jotter/directory/Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits', header=True)

fig3 = plt.figure()
plt.errorbar(x_340GHz - x_B3, y_340GHz - y_B3, xerr=data['x_err_340GHz'], yerr=data['y_err_340GHz'], c='b', alpha=0.4, label='340GHz', linestyle='', marker='*')
plt.errorbar(x_470GHz - x_B3, y_470GHz - y_B3, xerr=data['x_err_470GHz'], yerr=data['y_err_470GHz'], c='g', alpha=0.4, label='470GHz', linestyle='', marker='*')
plt.errorbar(x_B6 - x_B3, y_B6 - y_B3, xerr=data['x_err_B6'], yerr=data['y_err_B6'], c='y', alpha=0.4, label='B6', linestyle='', marker='*')
plt.scatter(RA_HC - x_B3, DEC_HC - y_B3, c='violet', alpha=0.4, label='HC IR')
plt.scatter(RA_RRS - x_B3, DEC_RRS - y_B3, c='k', alpha=0.4, label='RRS IR')

plt.xlabel('Delta RA (deg)')
plt.ylabel('Delta DEC (deg)')
plt.legend()
plt.grid()


'''
#Quiver plot of proper motions
mywcs = WCS(header).celestial
#img = img[::-1]
plt.imshow(img, vmin=0.5e-4, vmax=1e-3, origin='lower')

ind = np.where(np.isnan(data['RA_HC'])==False)
arrow_x = data['gauss_x_B3'][ind]
arrow_y = data['gauss_y_B3'][ind]

arrow_x_pix, arrow_y_pix = mywcs.all_world2pix(arrow_x, arrow_y, 1)

#arrow_x_pix = len(img) - arrow_x_pix
#arrow_y_pix = len(img) - arrow_y_pix

arrow_xlen = np.sin(((dist_table['pa_HC'][ind]-(np.pi/2))*u.degree).to(u.rad))*dist_table['v_HC'][ind]
arrow_ylen = -1*np.cos(((dist_table['pa_HC'][ind]-(np.pi/2))*u.degree).to(u.rad))*dist_table['v_HC'][ind]

for i in range(len(ind[0])):
	plt.text(arrow_x_pix[i], arrow_y_pix[i], data['_idx_B3'].data[ind][i], color='g')
plt.quiver(arrow_x_pix, arrow_y_pix, arrow_xlen, arrow_ylen, color='g')
'''
plt.show()


