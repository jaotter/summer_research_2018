import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, join, vstack
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from PM_fit import calc_pm
from calc_dates import *

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits'))

#MAX97 table
MAX97 = ascii.read('/lustre/aoc/students/jotter/tables/MAX97.txt')
MAX97['RA_MAX97'] = RA_to_deg(MAX97['RAh'], MAX97['RAm'], MAX97['RAs'])
MAX97['DEd'][np.where(MAX97['DE-'] == '-')] *= -1
MAX97['DEC_MAX97'] = DEC_to_deg(MAX97['DEd'], MAX97['DEm'], MAX97['DEs'])
MAX97.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
#MAX97['D_ID'] = np.full(len(MAX97), -1)
MAX97.rename_column('HC2000-Seq', 'HC_ID')
#MAX_coord = SkyCoord(MAX97['RA_MAX97']*u.degree, MAX97['DEC_MAX97']*u.degree)

MAX97['HC_ID'] = MAX97['HC_ID'].filled(-1).astype(int)
MAX_joined = join(data, MAX97, keys='HC_ID', join_type='left')

#tet Orion table
tet12 = ascii.read('/lustre/aoc/students/jotter/tables/tet2012.txt', data_start=4, header_start=2, guess=False, data_end=46, format='tab')
tet12['RA_tet12'] = RA_to_deg(5, 35, tet12['RA'])
tet_decmin = np.array([float(st[0:2]) for st in tet12['Decl']])
tet_decsec = np.array([float(st[3:-1]) for st in tet12['Decl']])
tet12['DEC_tet12'] = DEC_to_deg(-5, tet_decmin, tet_decsec)
tet12.remove_columns(('RA', 'Decl', 'X', 'Y'))
tet12['D_ID'] = np.full(len(tet12), -1)

tet_coord = SkyCoord(tet12['RA_tet12']*u.degree, tet12['DEC_tet12']*u.degree)
MAX_B3_coord = SkyCoord(MAX_joined['gauss_x_B3']*u.degree, MAX_joined['gauss_y_B3']*u.degree)

idx, d2d, d3d = tet_coord.match_to_catalog_sky(MAX_B3_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
tetmatch = np.where(d2d.degree < 0.5*(1/3600)) #matches are within 0.5 arcsec
for all_ind in tetmatch[0]:
	tet12[all_ind]['D_ID'] = MAX_joined[idx[all_ind]]['D_ID']

tetMAX_joined = join(MAX_joined, tet12, keys='D_ID', join_type='left')


#HH table
HH1 = ascii.read('/lustre/aoc/students/jotter/tables/HH508.txt', format='tab', data_start=1)
HH2 = ascii.read('/lustre/aoc/students/jotter/tables/HH508_2.txt', format='tab', data_start=1)
HH3 = ascii.read('/lustre/aoc/students/jotter/tables/HH508_3.txt', format='tab', data_start=1)
HH12 = vstack((HH1, HH2),join_type='outer')
HH = vstack((HH12, HH3),join_type='outer')
HH_rah =  np.array([int(st[0]) for st in HH['RA']])
HH_ram = np.array([int(st[2:4]) for st in HH['RA']])
HH_ras = np.array([float(st[5:-1]) for st in HH['RA']])
HH['RA_HH'] = RA_to_deg(HH_rah, HH_ram, HH_ras)
HH_ded =  np.array([int(st[0:2]) for st in HH['DEC']])
HH_dem = np.array([int(st[3:5]) for st in HH['DEC']])
HH_des = np.array([float(st[6:10]) for st in HH['DEC']])
HH['DEC_HH'] = DEC_to_deg(HH_ded, HH_dem, HH_des)
HH.remove_columns(('RA', 'DEC'))
HH['D_ID'] = np.full(len(HH), -1)

HH_coord = SkyCoord(HH['RA_HH']*u.degree, HH['DEC_HH']*u.degree)
MAX_tet_B3_coord = SkyCoord(tetMAX_joined['RA_HC']*u.degree, tetMAX_joined['DEC_HC']*u.degree)

idx, d2d, d3d = HH_coord.match_to_catalog_sky(MAX_tet_B3_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
HHmatch = np.where(d2d.degree < 1*(1/3600)) #matches are within 1 arcsec
for all_ind in HHmatch[0]:
	HH[all_ind]['D_ID'] = tetMAX_joined[idx[all_ind]]['D_ID']

HHtetMAX_joined = join(tetMAX_joined, HH, keys='D_ID', join_type='left')
'''
plt.scatter(HH['RA_HH'], HH['DEC_HH'], label='HH')
plt.scatter(MAX97['RA_MAX97'], MAX97['DEC_MAX97'], label='MAX97')
plt.scatter(tet12['RA_tet12'], tet12['DEC_tet12'], label='tet12')
plt.scatter(data['gauss_x_B3'], data['gauss_y_B3'], label='B3')
plt.scatter(data['RA_HC'], data['DEC_HC'], label='HC')
plt.legend()
plt.show()

HH508 = [5,35,16.05,-5,23,7.2]
MAX97 = [5,35,16.064,-5,23,7.13]
OriB3 = [5,35,16.064,-5,23,7.05]
OriB2 = [5,35,16.069,-5,23,6.96]
names = ['HH508_2000','MAX97_2004','OriB3_2012','OriB2_2012','HC2000','B3']
RAs = [RA_to_deg(HH508[0],HH508[1],HH508[2]), RA_to_deg(MAX97[0],MAX97[1],MAX97[2]),RA_to_deg(OriB3[0],OriB3[1],OriB3[2]),RA_to_deg(OriB2[0],OriB2[1],OriB2[2]),data['RA_HC'][60],data['gauss_x_B3'][60]]
DECs = [DEC_to_deg(HH508[3],HH508[4],HH508[5]), DEC_to_deg(MAX97[3],MAX97[4],MAX97[5]),DEC_to_deg(OriB3[3],OriB3[4],OriB3[5]),DEC_to_deg(OriB2[3],OriB2[4],OriB2[5]),data['DEC_HC'][60],data['gauss_y_B3'][60]]

for i in range(len(RAs)):
	plt.scatter(RAs[i], DECs[i], label=names[i])
plt.legend()

'''

ind0 = np.where(HHtetMAX_joined['_idx_B3'] == 0)[0][0] 
ind_ref = np.where(HHtetMAX_joined['_idx_B3'] == 3)

HC_ra_err = (0.017*u.arcsecond).to(u.degree)
HC_dec_err = (0.028*u.arcsecond).to(u.degree)

RAs = [float(tetMAX_joined['RA_tet12'][ind0] - tetMAX_joined['RA_tet12'][ind_ref]), float(tetMAX_joined['RA_MAX97'][ind0] - tetMAX_joined['RA_MAX97'][ind_ref]), float(tetMAX_joined['RA_HC'][ind0] - tetMAX_joined['RA_HC'][ind_ref]), float(tetMAX_joined['gauss_x_B3'][ind0] - tetMAX_joined['gauss_x_B3'][ind_ref]), float(tetMAX_joined['RA_Fb'][ind0] - tetMAX_joined['RA_Fb'][ind_ref])]
DECs = [float(tetMAX_joined['DEC_tet12'][ind0] - tetMAX_joined['DEC_tet12'][ind_ref]), float(tetMAX_joined['DEC_MAX97'][ind0] - tetMAX_joined['DEC_MAX97'][ind_ref]), float(tetMAX_joined['DEC_HC'][ind0] - tetMAX_joined['DEC_HC'][ind_ref]), float(tetMAX_joined['gauss_y_B3'][ind0] - tetMAX_joined['gauss_y_B3'][ind_ref]), float(tetMAX_joined['DEC_Fb'][ind0] - tetMAX_joined['DEC_Fb'][ind_ref])]
names = ['tet2012 (2011)','MAX97 (2002)','HC2000 (2000)','B3 (2018)','Forbrich (2012)']
RA_err = [0.1/(3600),0.1/(3600),HC_ra_err.value,tetMAX_joined['x_err_B3'][ind0],0.006/(3600)]
DEC_err = [0.1/(3600),0.1/(3600),HC_dec_err.value,tetMAX_joined['y_err_B3'][ind0],0.007/(3600)]

times = [12,2+1/12,0,18+1/12,13]
times_err = [0.001,1/12,0.001,0.001,0.01]

v = calc_pm(RAs, DECs, RA_err, DEC_err, times, times_err)
v = v*u.degree/u.year
print(v.to(u.mas/u.year))
plt.show()
'''
for i in range(len(RAs)):
	plt.errorbar(RAs[i], DECs[i], xerr=RA_err[i], yerr=DEC_err[i], label=names[i], linestyle='')

plt.xlabel('RA (deg)')
plt.ylabel('DEC (deg)')
plt.grid()
plt.legend()


def f(m, x): #m is list of parameters (slope and intercept)
	return m[0]*x + m[1]

linear = Model(f)
mydata = RealData(RAs, DECs, sx=RA_err, sy=DEC_err)
myodr = ODR(mydata, linear, beta0=[1,0])
myoutput=myodr.run()

x_vals = np.linspace(np.min(RAs), np.max(RAs),10)
plt.plot(x_vals, f(myoutput.beta, x_vals))

plt.show() '''
