from astropy.io import fits, ascii
from astropy.table import Table, join
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np

def RA_to_deg(HH, MM, SS):
	return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
	return DD + (MM/60 + SS/(60*60))*np.sign(DD)

def RA_to_hr(deg):
    hr = deg*(24/360)
    minute = (hr - int(hr))*60
    sec = (minute - int(minute))*60
    return int(hr), int(minute), sec

def DEC_to_arcsec(deg):
    degree = int(deg)
    minute = np.sign(deg)*(deg - int(deg))*60
    sec = (minute - int(minute))*60
    return degree, int(minute), sec


B3data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B3_500klplus_allsrcs_catalog_B3_ref.txt')
B6data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B6_500klplus_allsrcs_catalog_B3_ref.txt')
B7data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B7_500klplus_allsrcs_catalog_B3_ref.txt')
A340data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/A340_500klplus_allsrcs_catalog_B3_ref.txt')
A470data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/A470_500klplus_allsrcs_catalog_B3_ref.txt')


B3B6 = join(B3data, B6data, keys='D_ID', join_type = 'outer')
B3B6B7 = join(B3B6, B7data, keys='D_ID', join_type = 'outer')
B3B6B7A340 = join(B3B6B7, A340data, keys='D_ID', join_type = 'outer')
all_bands = join(B3B6B7A340, A470data, keys='D_ID', join_type = 'outer')

all_bands_coord = SkyCoord(ra=all_bands['RA_B3']*u.degree, dec=all_bands['DEC_B3']*u.degree)

#loading in the Forbrich catalog
Fb_data = ascii.read('/lustre/aoc/students/jotter/tables/Forbrich2016.txt')
Fb_data.rename_column('Seq', 'Fb_ID')
Fb_data.rename_column('e_RAs', 'e_RAs_Fb')
Fb_data.rename_column('e_DEs', 'e_DEs_Fb')
Fb_data.rename_column('Speak', 'Speak_Fb')
Fb_data.rename_column('e_Speak', 'e_Speak_Fb')
#converting from hours,minutes,seconds to degrees
ra_Fb = RA_to_deg(Fb_data['RAh'],Fb_data['RAm'],Fb_data['RAs'])*u.degree
#degrees + arcmins and arcsecs
Fb_data['DEd'][np.where(Fb_data['DE-'] == '-')] *= -1
dec_Fb = DEC_to_deg(Fb_data['DEd'],Fb_data['DEm'],Fb_data['DEs'])*u.degree
Fb_data['DEC_Fb'] = dec_Fb
Fb_data['RA_Fb'] = ra_Fb
Fb_data.remove_columns(['RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs'])
ID_empty_Fb = np.full(len(Fb_data),-1)
Fb_data['D_ID'] = ID_empty_Fb

Fb_coord = SkyCoord(ra=Fb_data['RA_Fb'], dec=Fb_data['DEC_Fb'])

#match Forbrich with dendrogram data
idx, d2d, d3d = all_bands_coord.match_to_catalog_sky(Fb_coord)
#only qualify as a match if they are within 0.5"
matches = np.where(d2d.degree < 0.5*(1/3600))

for all_ind in matches[0]:
	Fb_data[idx[all_ind]]['D_ID'] = all_bands[all_ind]['D_ID']

Fb_dend_table = join(all_bands, Fb_data, keys='D_ID', join_type='left')
Fb_dend_coord = SkyCoord(ra=Fb_dend_table['RA_B3']*u.degree, dec=Fb_dend_table['DEC_B3']*u.degree)


#loading in the HC catalog
HC_data = ascii.read('/lustre/aoc/students/jotter/tables/HC2000.txt')
HC_data.rename_column('ID', 'HC_ID')
ID_empty_HC = np.full(len(HC_data),'N/A',dtype='U10') #add column for RRS ID to match later
HC_data['RRS_ID'] = ID_empty_HC
#converting from hours,minutes,seconds to degrees
ra_HC = RA_to_deg(HC_data['RAh'],HC_data['RAm'],HC_data['RAs'])*u.degree
#degrees + arcmins and arcsecs
HC_data['DEd'][np.where(HC_data['DE-'] == '-')] *= -1
dec_HC = DEC_to_deg(HC_data['DEd'],HC_data['DEm'],HC_data['DEs'])*u.degree
HC_data['DEC_HC'] = dec_HC
HC_data['RA_HC'] = ra_HC
HC_data['D_ID'] = np.full(len(HC_data), -1)
HC_data.remove_columns(['RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs'])

HC_table = Table([HC_data['HC_ID'],HC_data['RA_HC'],HC_data['DEC_HC'],HC_data['D_ID']])
HC_coord = SkyCoord(ra=HC_table['RA_HC'], dec=HC_table['DEC_HC'])


idx, d2d, d3d = Fb_dend_coord.match_to_catalog_sky(HC_coord)
#only qualify as a match if they are within 0.5"
matches = np.where(d2d.degree < 0.5*(1/3600))

for all_ind in matches[0]:
	HC_table[idx[all_ind]]['D_ID'] = Fb_dend_table[all_ind]['D_ID']

full_table = join(Fb_dend_table, HC_table, keys='D_ID', join_type='left')





