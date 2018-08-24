from astropy.io import fits, ascii
from astropy.table import Table, join, vstack
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from PM_fit import calc_pm
from calc_dates import *
#currently 11 catalogs
COUP_poserr = (np.nan,np.nan) #these have source-dependent errors, so deal with them later
FDM_poserr = (np.nan,np.nan)
MAX_poserr = (np.nan,np.nan)
B3_poserr = (np.nan,np.nan)
B6_poserr = (np.nan,np.nan)
B7_poserr = (np.nan,np.nan)
A340_poserr = (np.nan,np.nan)
A470_poserr = (np.nan,np.nan)

HC_poserr = (0.017/3600, 0.028/3600)
Fb_poserr = (0.006/3600, 0.007/3600)
MLLA_poserr = (0.1/3600, 0.1/3600) #MLLA error about 0.1"
RRS_poserr = (0.025/3600, 0.025/3600) #25mas, mentioned in paper but not certain
MRD_poserr = (0.1/3600, 0.1/3600) #can't find actual error, conservative estimate
OW94_poserr = (0.1/3600, 0.1/3600) #can't find actual error, conservative estimate
FW_poserr = (0.1/3600, 0.1/3600)
tet12_poserr = (0.1/3600, 0.1/3600)
LML_poserr = (0.1/3600, 0.1/3600)
HH_poserr = (0.1/3600, 0.1/3600)
GBS_poserr = (0.1/3600, 0.1/3600)
MGM_poserr = (0.1/3600, 0.1/3600)

#Observation specifications of various catalogs
obs = Table(names=('name', 'month', 'year', 't_err', 'ra_err', 'dec_err'), dtype=('S6', 'i2', 'i4', 'f8', 'f8', 'f8'))
FW_date, FW_err = middle_date(12,2007,4,2008)
obs.add_row(('FW', FW_date[0], FW_date[1], FW_err, FW_poserr[0], FW_poserr[1]))
HC_date = (10,1999)
HC_err = 0.5/12
obs.add_row(('HC', HC_date[0], HC_date[1], HC_err, HC_poserr[0], HC_poserr[1]))
B3_date = (10,2017)
B3_err = 0.5/12
obs.add_row(('B3', B3_date[0], B3_date[1], B3_err, B3_poserr[0], B3_poserr[1]))
B6_date = (10,2017)
B6_err = 0.5/12
obs.add_row(('B6', B6_date[0], B6_date[1], B6_err, B6_poserr[0], B6_poserr[1]))
B7_date = (10,2017)
B7_err = 0.5/12
obs.add_row(('B7', B7_date[0], B7_date[1], B7_err, B7_poserr[0], B7_poserr[1]))
A340_date = (11,2017)
A340_err = 0.5/12
obs.add_row(('A340', A340_date[0], A340_date[1], A340_err, A340_poserr[0], A340_poserr[1]))
A470_date = (8,2015)
A470_err = 0.5/12
obs.add_row(('A470', A470_date[0], A470_date[1], A470_err, A470_poserr[0], A470_poserr[1]))
LML_date, LML_err = middle_date(10,2002,12,2002)
obs.add_row(('LML', LML_date[0], LML_date[1], LML_err, LML_poserr[0], LML_poserr[1]))
MLLA_date, MLLA_err = middle_date(12,1997,3,2000)
obs.add_row(('MLLA', MLLA_date[0], MLLA_date[1], MLLA_err, MLLA_poserr[0], MLLA_poserr[1]))
OW94_date = (1,1994)
OW94_err = 0.5/12
obs.add_row(('OW94', OW94_date[0], OW94_date[1], OW94_err, OW94_poserr[0], OW94_poserr[1]))
COUP_date = (1,2003)
COUP_err = 0.5/12
obs.add_row(('COUP', COUP_date[0], COUP_date[1], COUP_err, COUP_poserr[0], COUP_poserr[1]))
FDM_date = (2,2000)
FDM_err = 0.5/12
obs.add_row(('FDM', FDM_date[0], FDM_date[1], FDM_err, FDM_poserr[0], FDM_poserr[1]))
tet12_date = (10,2011)
tet12_err = 0.5/12
obs.add_row(('tet12', tet12_date[0], tet12_date[1], tet12_err, tet12_poserr[0], tet12_poserr[1]))
Fb_date = (10,2012)
Fb_err = 0.5/12
obs.add_row(('Fb', Fb_date[0], Fb_date[1], Fb_err, Fb_poserr[0], Fb_poserr[1]))
HH_date, HH_err = middle_date(4,1998,10,1999)
obs.add_row(('HH', HH_date[0], HH_date[1], HH_err, HH_poserr[0], HH_poserr[1]))
MAX_date, MAX_err = middle_date(11,1998,12,2000)
obs.add_row(('MAX', MAX_date[0], MAX_date[1], MAX_err, MAX_poserr[0], MAX_poserr[1]))
MRD_date, MRD_err = middle_date(10,2004,4,2005)
obs.add_row(('MRD', MRD_date[0], MRD_date[1], MRD_err, MRD_poserr[0], MRD_poserr[1]))
RRS_date, RRS_err = middle_date(11,2004,4,2005)
obs.add_row(('RRS', RRS_date[0], RRS_date[1], RRS_err, RRS_poserr[0], RRS_poserr[1]))
MGM_date = (10,2004)
MGM_err = 0.5/12
obs.add_row(('MGM', MGM_date[0], MGM_date[1], MGM_err, MGM_poserr[0], MGM_poserr[1]))
GBS_date, GBS_err = middle_date(7,2011,9,2011)
obs.add_row(('GBS', GBS_date[0], GBS_date[1], GBS_err, GBS_poserr[0], GBS_poserr[1]))

obs.write('/lustre/aoc/students/jotter/dendro_catalogs/obs_dates_errs.txt', format='ascii', overwrite=True)

#Matching all the misc catalogs for PM measurements

data = Table(fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits')) 

#MAX table
MAX = ascii.read('/lustre/aoc/students/jotter/tables/MAX.txt')
MAX['RA_MAX'] = RA_to_deg(MAX['RAh'], MAX['RAm'], MAX['RAs'])
MAX['DEd'][np.where(MAX['DE-'] == '-')] *= -1
MAX['DEC_MAX'] = DEC_to_deg(MAX['DEd'], MAX['DEm'], MAX['DEs'])
MAX.rename_column('HC2000-Seq', 'HC_ID')
MAX.rename_column('Seq', 'MAX_ID')
MAX.rename_column('PosFlag', 'PosFlag_MAX')
MAX['HC_ID'] = MAX['HC_ID'].filled(-1).astype(int)
MAX_table = Table((MAX['RA_MAX'], MAX['DEC_MAX'], MAX['HC_ID'], MAX['MAX_ID'], MAX['PosFlag_MAX']))
data = join(data, MAX_table, keys='HC_ID', join_type='left')

#MRD2012 data
MRD = ascii.read('/lustre/aoc/students/jotter/tables/MRD2012-2.txt')
MRD['RA_MRD'] = RA_to_deg(MRD['RAh'], MRD['RAm'], MRD['RAs'])
MRD['DEd'][np.where(MRD['DE-'] == '-')] *= -1
MRD['DEC_MRD'] = DEC_to_deg(MRD['DEd'], MRD['DEm'], MRD['DEs'])
MRD.rename_column('OM', 'MRD_ID')
MRD['D_ID'] = np.full(len(MRD), -1)

MRD_table = Table((MRD['RA_MRD'], MRD['DEC_MRD'], MRD['D_ID'], MRD['MRD_ID']))

MRD_coord = SkyCoord(MRD_table['RA_MRD'].quantity.value*u.degree, MRD_table['DEC_MRD'].quantity.value*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)
idx, d2d, d3d = MRD_coord.match_to_catalog_sky(HC_coord)
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec
for all_ind in match[0]:
	MRD_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, MRD_table, keys='D_ID', join_type='left')

#MGM2012 data
MGM = ascii.read('/lustre/aoc/students/jotter/tables/MGM2012.txt')
MGM['RA_MGM'] = RA_to_deg(MGM['RAh'], MGM['RAm'], MGM['RAs'])
MGM['DEd'][np.where(MGM['DE-'] == '-')] *= -1
MGM['DEC_MGM'] = DEC_to_deg(MGM['DEd'], MGM['DEm'], MGM['DEs'])
MGM.rename_column('Num', 'MGM_ID')
MGM['D_ID'] = np.full(len(MGM), -1)

MGM_table = Table((MGM['RA_MGM'], MGM['DEC_MGM'], MGM['D_ID'], MGM['MGM_ID']))

MGM_coord = SkyCoord(MGM_table['RA_MGM'].quantity.value*u.degree, MGM_table['DEC_MGM'].quantity.value*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)
idx, d2d, d3d = MGM_coord.match_to_catalog_sky(HC_coord)
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec
for all_ind in match[0]:
	MGM_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, MGM_table, keys='D_ID', join_type='left')


#OW94 table
OW94 = ascii.read('/lustre/aoc/students/jotter/tables/OW94.txt', data_start=7, data_end=387, header_start=4)
OW94['RA_OW94'] = RA_to_deg(OW94['RAh'], OW94['RAm'], OW94['RAs'])
OW94['DEC_OW94'] = DEC_to_deg(OW94['DEd'], OW94['DEm'], OW94['DEs'])
OW94.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
OW94.rename_column('[OW94]', 'OW94_ID')
OW94['D_ID'] = np.full(len(OW94), -1)

OW_coord = SkyCoord(OW94['RA_OW94']*u.degree, OW94['DEC_OW94']*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)

idx, d2d, d3d = OW_coord.match_to_catalog_sky(HC_coord)
match = np.where(d2d.degree < 2*(1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	OW94[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, OW94, keys='D_ID', join_type='left')

#MLLA table
MLLA = ascii.read('/lustre/aoc/students/jotter/tables/MLLA.txt')
MLLA['RA_MLLA'] = RA_to_deg(MLLA['RAh'], MLLA['RAm'], MLLA['RAs'])
MLLA['DEd'][np.where(MLLA['DE-'] == '-')] *= -1
MLLA['DEC_MLLA'] = DEC_to_deg(MLLA['DEd'], MLLA['DEm'], MLLA['DEs'])
MLLA.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
MLLA.rename_column('HC2000-Seq', 'HC_ID')
MLLA.rename_column('Seq', 'MLLA_ID')
MLLA['HC_ID'] = MLLA['HC_ID'].filled(-1).astype(int)
MLLA_table = Table((MLLA['RA_MLLA'], MLLA['DEC_MLLA'], MLLA['MLLA_ID'], MLLA['HC_ID']))
data = join(data, MLLA_table, keys='HC_ID', join_type='left')

#COUP data
COUP = ascii.read('/lustre/aoc/students/jotter/tables/COUP.txt')
COUP.rename_column('RAdeg', 'RA_COUP')
COUP.rename_column('DEdeg', 'DEC_COUP')
COUP.rename_column('PosErr', 'PosErr_COUP')
COUP.rename_column('Seq', 'COUP_ID')
COUP['D_ID'] = np.full(len(COUP), -1)

COUP_table = Table((COUP['RA_COUP'], COUP['DEC_COUP'], COUP['PosErr_COUP'], COUP['COUP_ID'], COUP['D_ID']))

COUP_coord = SkyCoord(COUP_table['RA_COUP'], COUP_table['DEC_COUP'])
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)

idx, d2d, d3d = COUP_coord.match_to_catalog_sky(HC_coord)
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	COUP_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, COUP_table, keys='D_ID', join_type='left')

#LML table
LML = ascii.read('/lustre/aoc/students/jotter/tables/LML2004.txt')
LML['RA_LML'] = RA_to_deg(LML['RAh'], LML['RAm'], LML['RAs'])
LML['DEd'][np.where(LML['DE-'] == '-')] *= -1
LML['DEC_LML'] = DEC_to_deg(LML['DEd'], LML['DEm'], LML['DEs'])
LML.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
LML.rename_column('HC2000-Seq', 'HC_ID')
LML.rename_column('Seq', 'LML_ID')
LML.rename_column('PosFlag', 'PosFlag_LML')
LML['HC_ID'] = LML['HC_ID'].filled(-1).astype(int)
LML_table = Table((LML['RA_LML'], LML['DEC_LML'], LML['HC_ID'], LML['LML_ID'], LML['PosFlag_LML']))
data = join(data, LML_table, keys='HC_ID', join_type='left')

#FW2011 table
FW2011 = ascii.read('/lustre/aoc/students/jotter/tables/FW2011.txt', data_start=3, data_end=34, header_start=2)
RAh = np.array([int(st[0:2]) for st in FW2011['alpha_J2000']])
RAm = np.array([int(st[3:5]) for st in FW2011['alpha_J2000']])
RAs = np.array([float(st[6:12]) for st in FW2011['alpha_J2000']])
FW2011['RA_FW'] = RA_to_deg(RAh, RAm, RAs)
DEh = np.array([int(st[0:3]) for st in FW2011['delta_J2000']])
DEm = np.array([int(st[4:6]) for st in FW2011['delta_J2000']])
DEs = np.array([float(st[7:13]) for st in FW2011['delta_J2000']])
FW2011['DEC_FW'] = DEC_to_deg(DEh, DEm, DEs)
FW2011.remove_columns(('alpha_J2000', 'delta_J2000'))
FW2011['D_ID'] = np.full(len(FW2011), -1)
FW2011.rename_column('ID', 'ID_FW')

FW_coord = SkyCoord(FW2011['RA_FW']*u.degree, FW2011['DEC_FW']*u.degree)
B3_coord = SkyCoord(data['gauss_x_B3']*u.degree, data['gauss_y_B3']*u.degree)

idx, d2d, d3d = FW_coord.match_to_catalog_sky(B3_coord)
match = np.where(d2d.degree < 0.5*(1/3600)) #matches are within 1 arcsec
for all_ind in match[0]:
	FW2011[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']
data = join(data, FW2011, keys='D_ID', join_type='left')

#tet Orion table
tet12 = ascii.read('/lustre/aoc/students/jotter/tables/tet2012.txt', data_start=4, header_start=2, guess=False, data_end=46, format='tab')
tet12['RA_tet12'] = RA_to_deg(5, 35, tet12['RA'])
tet_decmin = np.array([float(st[0:2]) for st in tet12['Decl']])
tet_decsec = np.array([float(st[3:-1]) for st in tet12['Decl']])
tet12['DEC_tet12'] = DEC_to_deg(-5, tet_decmin, tet_decsec)
tet12['D_ID'] = np.full(len(tet12), -1)

tet12_table = Table((tet12['RA_tet12'], tet12['DEC_tet12'], tet12['D_ID']))

tet_coord = SkyCoord(tet12_table['RA_tet12']*u.degree, tet12_table['DEC_tet12']*u.degree)
B3_coord = SkyCoord(data['gauss_x_B3']*u.degree, data['gauss_y_B3']*u.degree)

idx, d2d, d3d = tet_coord.match_to_catalog_sky(B3_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
tetmatch = np.where(d2d.degree < 0.5*(1/3600)) #matches are within 0.5 arcsec
for all_ind in tetmatch[0]:
	tet12_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, tet12_table, keys='D_ID', join_type='left')


#GBS table
GBS = ascii.read('/lustre/aoc/students/jotter/tables/GBS.txt')
GBS_rah =  np.array([int(st[1:3]) for st in GBS['ID']])
GBS_ram = np.array([int(st[3:5]) for st in GBS['ID']])
GBS_ras = np.array([float(st[5:10]) for st in GBS['ID']])
GBS['RA_GBS'] = RA_to_deg(GBS_rah, GBS_ram, GBS_ras)
GBS_ded =  np.array([int(st[10:13]) for st in GBS['ID']])
GBS_dem = np.array([int(st[13:15]) for st in GBS['ID']])
GBS_des = np.array([float(st[15:19]) for st in GBS['ID']])
GBS['DEC_GBS'] = DEC_to_deg(GBS_ded, GBS_dem, GBS_des)
GBS['D_ID'] = np.full(len(GBS), -1)
GBS.rename_column('ID', 'GBS_ID')
GBS_table = Table((GBS['RA_GBS'], GBS['DEC_GBS'], GBS['D_ID'], GBS['GBS_ID']))

GBS_coord = SkyCoord(GBS_table['RA_GBS']*u.degree, GBS_table['DEC_GBS']*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)

idx, d2d, d3d = GBS_coord.match_to_catalog_sky(HC_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
GBSmatch = np.where(d2d.degree < 5*(1/3600)) #matches are within 1 arcsec
for all_ind in GBSmatch[0]:
	GBS_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, GBS_table, keys='D_ID', join_type='left')


#HH table
HH1 = ascii.read('/lustre/aoc/students/jotter/tables/HH508.txt', format='tab', data_start=1)
HH2 = ascii.read('/lustre/aoc/students/jotter/tables/HH508_2.txt', format='tab', data_start=1)
HH3 = ascii.read('/lustre/aoc/students/jotter/tables/HH508_3.txt', format='tab', data_start=1)
HH12 = vstack((HH1, HH2),join_type='outer')
HH = vstack((HH12, HH3),join_type='outer')
HH_rah =  np.array([int(st[0]) for st in HH['RA']])
HH_ram = np.array([int(st[2:4]) for st in HH['RA']])
HH_ras = np.array([float(st[5:10]) for st in HH['RA']])
HH['RA_HH'] = RA_to_deg(HH_rah, HH_ram, HH_ras)
HH_ded =  np.array([int(st[0:2]) for st in HH['DEC']])
HH_dem = np.array([int(st[3:5]) for st in HH['DEC']])
HH_des = np.array([float(st[6:10]) for st in HH['DEC']])
HH['DEC_HH'] = DEC_to_deg(HH_ded, HH_dem, HH_des)
HH['D_ID'] = np.full(len(HH), -1)
HH.rename_column('HH No.', 'HH_ID')
HH_table = Table((HH['RA_HH'], HH['DEC_HH'], HH['D_ID'], HH['HH_ID']))

HH_coord = SkyCoord(HH_table['RA_HH']*u.degree, HH_table['DEC_HH']*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)

idx, d2d, d3d = HH_coord.match_to_catalog_sky(HC_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
HHmatch = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec
for all_ind in HHmatch[0]:
	HH_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, HH_table, keys='D_ID', join_type='left')

#FDM table
FDM = ascii.read('/lustre/aoc/students/jotter/tables/FDM2003.txt')
FDM['RA_FDM'] = RA_to_deg(FDM['RAh'], FDM['RAm'], FDM['RAs'])
FDM['DEd'][np.where(FDM['DE-'] == '-')] *= -1
FDM['DEC_FDM'] = DEC_to_deg(FDM['DEd'], FDM['DEm'], FDM['DEs'])
FDM.remove_columns(('RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs'))
FDM['D_ID'] = np.full(len(FDM), -1)
FDM.rename_column('PosErr', 'PosErr_FDM')
FDM.rename_column('Num', 'FDM_ID')
FDM_table = Table((FDM['RA_FDM'], FDM['DEC_FDM'], FDM['FDM_ID'], FDM['PosErr_FDM'], FDM['D_ID']))

FDM_coord = SkyCoord(FDM_table['RA_FDM'].quantity.value*u.degree, FDM_table['DEC_FDM'].quantity.value*u.degree)
HC_coord = SkyCoord(data['RA_HC']*u.degree, data['DEC_HC']*u.degree)

idx, d2d, d3d = FDM_coord.match_to_catalog_sky(HC_coord)
#idx is a list of indices for data with the list index corresponding to the match in data
match = np.where(d2d.degree < (1/3600)) #matches are within 1 arcsec

for all_ind in match[0]:
	FDM_table[all_ind]['D_ID'] = data[idx[all_ind]]['D_ID']

data = join(data, FDM_table, keys='D_ID', join_type='left')

data.rename_column('gauss_x_B3', 'RA_B3')
data.rename_column('gauss_y_B3', 'DEC_B3')
data.rename_column('gauss_x_B6', 'RA_B6')
data.rename_column('gauss_y_B6', 'DEC_B6')
#data.rename_column('gauss_x_B7', 'RA_B7')
#data.rename_column('gauss_y_B7', 'DEC_B7')
data.rename_column('gauss_x_340GHz', 'RA_A340')
data.rename_column('gauss_y_340GHz', 'DEC_A340')
data.rename_column('gauss_x_470GHz', 'RA_A470')
data.rename_column('gauss_y_470GHz', 'DEC_A470')

data.rename_column('x_err_340GHz', 'x_err_A470')
data.rename_column('y_err_340GHz', 'y_err_A470')

data.write('/lustre/aoc/students/jotter/dendro_catalogs/simbad_catalog.fits',format='fits',overwrite=True)



