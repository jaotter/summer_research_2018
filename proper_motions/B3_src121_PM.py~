import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from PM_fit import calc_pm

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

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'))

OW94 = ascii.read('/lustre/aoc/students/jotter/tables/OW94.txt', data_start=7, data_end=387, header_start=4)
OW94['RA_OW94'] = RA_to_deg(OW94['RAh'], OW94['RAm'], OW94['RAs'])
OW94['DEC_OW94'] = DEC_to_deg(OW94['DEd'], OW94['DEm'], OW94['DEs'])
OW94.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
OW94['D_ID'] = np.full(len(OW94), -1)

OW_coord = SkyCoord(OW94['RA_OW94']*u.degree, OW94['DEC_OW94']*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)


idx, d2d, d3d = OW_coord.match_to_catalog_sky(HC_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	OW94[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

OW_joined = join(data, OW94, keys='D_ID', join_type='left')

#MRD2012 data
MRD = ascii.read('/lustre/aoc/students/jotter/tables/MRD2012.txt')
MRD['RA_MRD'] = RA_to_deg(MRD['RAh'], MRD['RAm'], MRD['RAs'])
MRD['DEd'][np.where(MRD['DE-'] == '-')] *= -1
MRD['DEC_MRD'] = DEC_to_deg(MRD['DEd'], MRD['DEm'], MRD['DEs'])
#MRD.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
MRD['D_ID'] = np.full(len(MRD), -1)
MRD_coord = SkyCoord(MRD['RA_MRD'].quantity.value*u.degree, MRD['DEC_MRD'].quantity.value*u.degree)

HC_coord = SkyCoord(OW_joined['RA_HC']*u.degree, OW_joined['DEC_HC']*u.degree)
idx, d2d, d3d = MRD_coord.match_to_catalog_sky(HC_coord)
match = np.where(d2d.degree < 2*(1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	MRD[all_ind]['D_ID'] = OW_joined[idx[all_ind]]['D_ID']

MRD_OW = join(OW_joined, MRD, keys='D_ID', join_type='left')

#COUP data
COUP = ascii.read('/lustre/aoc/students/jotter/tables/COUP.txt')
COUP.rename_column('RAdeg', 'RA_COUP')
COUP.rename_column('DEdeg', 'DEC_COUP')
COUP['D_ID'] = np.full(len(COUP), -1)
COUP_coord = SkyCoord(COUP['RA_COUP'].quantity.value*u.degree, COUP['DEC_COUP'].quantity.value*u.degree)

HC_coord = SkyCoord(OW_joined['RA_HC']*u.degree, OW_joined['DEC_HC']*u.degree)
idx, d2d, d3d = COUP_coord.match_to_catalog_sky(HC_coord)
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	COUP[all_ind]['D_ID'] = MRD_OW[idx[all_ind]]['D_ID']

COUP_MRD_OW = join(MRD_OW, COUP, keys='D_ID', join_type='left')


#plt.scatter(OW94['RA_OW94'], OW94['DEC_OW94'], label='OW94')
#plt.scatter(MRD['RA_MRD'], MRD['DEC_MRD'], label='MRD2012')
#plt.scatter(MRD_OW['RA_HC'], MRD_OW['DEC_HC'], label='HC2000')
#plt.legend()
#plt.show()

ind83 = np.where(MRD_OW['_idx_B3'] == 83)
ind_ref = np.where(MRD_OW['_idx_B3'] == 0)

HC_ra_err = (0.017*u.arcsecond).to(u.degree)
HC_dec_err = (0.028*u.arcsecond).to(u.degree)

RAs = [float(MRD_OW['RA_OW94'][ind83] - MRD_OW['RA_OW94'][ind_ref]),float(MRD_OW['RA_HC'][ind83] - MRD_OW['RA_HC'][ind_ref]),float(MRD_OW['gauss_x_B3'][ind83] - MRD_OW['gauss_x_B3'][ind_ref]), float(MRD_OW['RA_MRD'][ind83] - MRD_OW['RA_MRD'][ind_ref])]
DECs = [float(MRD_OW['DEC_OW94'][ind83] - MRD_OW['DEC_OW94'][ind_ref]),float(MRD_OW['DEC_HC'][ind83] - MRD_OW['DEC_HC'][ind_ref]),float(MRD_OW['gauss_y_B3'][ind83] - MRD_OW['gauss_y_B3'][ind_ref]),float(MRD_OW['DEC_MRD'][ind83] - MRD_OW['DEC_MRD'][ind_ref])]
RA_err = [0.1/(3600),HC_ra_err.value,MRD_OW['x_err_B3'][ind83], 0.1/3600]
DEC_err = [0.1/(3600),HC_dec_err.value,MRD_OW['y_err_B3'][ind83], 0.1/3600]
times=[0,5+9/12,23+10/12,11]
times_err=[0.01,0.002, 0.001,3/12]

v = calc_pm(RAs, DECs, RA_err, DEC_err, times, times_err)
plt.show()
