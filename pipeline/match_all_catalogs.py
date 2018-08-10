#This file matches the dendrogram catalogs with the HC2000 catalog and the RRS2008 catalog
from astropy.io import fits, ascii
from astropy import units as u
import numpy as np
from astropy.coordinates import SkyCoord, Angle
import matplotlib.pyplot as plt
from astropy.table import Table,hstack,vstack,join,Column
import csv
import glob

def RA_to_deg(HH, MM, SS):
	return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
	return DD + (MM/60 + SS/(60*60))*np.sign(DD)

#first need to join together dendrogram catalogs
dendro_cats = sorted(glob.glob('/lustre/aoc/students/jotter/dendro_catalogs/*_dendro_catalog_leaves_fixed.fits'))
names = ['340GHz', '470GHz', 'B3', 'B6', 'B7_hr', 'B7_lr']

#dendro_table is the master table, start by setting it to the first table and add some columns
dendro_table = Table(fits.getdata(dendro_cats[0]))
dendro_table['D_ID'] = np.arange(0,len(dendro_table)) #column of 'master IDs', common to all dendrogram catalogs
dendro_table['RA'] = dendro_table['gauss_x_'+names[0]] #start master RA/Dec columns
dendro_table['DEC'] = dendro_table['gauss_y_'+names[0]]
dendro_table.remove_columns(['area_ellipse', 'area_exact', 'major_sigma', 'minor_sigma', 'position_angle', 'radius', 'x_cen', 'y_cen'])

#loop through other catalogs and add their data to dendro_table
for ind in np.arange(1,len(dendro_cats)):
	print(names[ind])
	#need to match against all previous sources
	ra_all = dendro_table['RA']
	dec_all = dendro_table['DEC']
	coord_all = SkyCoord(ra=ra_all*u.degree, dec=dec_all*u.degree)

	cat = Table(fits.getdata(dendro_cats[ind])) #cat is the next catalog to add to dendro_table
	cat['D_ID'] = np.full(len(cat),-1)
	if names[ind] != 'B7_hr' and names[ind] != 'B7_lr':
		cat.remove_columns(['area_ellipse', 'area_exact', 'major_sigma', 'minor_sigma', 'position_angle', 'radius', 'x_cen', 'y_cen'])
	ra2 = cat['gauss_x_'+names[ind]]
	dec2 = cat['gauss_y_'+names[ind]]
	coord2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
	idx, d2d, d3d = coord_all.match_to_catalog_sky(coord2)
	#idx is a list of indices for cat with the list index corresponding to the match in dendro_table
	match = np.where(d2d.degree < 0.1*(1/3600)) #matches are within 0.1 arcsec
	#print(names[ind])
	#match is indices of idx of matches
	#now loop through all sources and change the D_ID
	for all_ind in match[0]:
		cat[idx[all_ind]]['D_ID'] = dendro_table[all_ind]['D_ID']
	
		
	dendro_table = join(dendro_table, cat, keys='D_ID', join_type='outer', table_names=[names[ind -1], names[ind]])
	#print(dendro_table[match]['_idx_'+names[ind-1]])
	#give unmatched sources a new, unique D_ID
	unmatched = np.where(dendro_table['D_ID'] == -1)
	for un in unmatched[0]:
		dendro_table[un]['D_ID'] = np.nanmax(np.concatenate((dendro_table['D_ID'], cat['D_ID'])))+1
		cat_ind = np.where(cat['_idx_'+names[ind]] == dendro_table['_idx_'+names[ind]][un])
		dendro_table[un]['RA'] = cat[cat_ind]['gauss_x_' + names[ind]]
		dendro_table[un]['DEC'] = cat[cat_ind]['gauss_y_' + names[ind]]

dendro_coord = 	SkyCoord(ra=dendro_table['RA']*u.degree, dec=dendro_table['DEC']*u.degree)

dendro_table.write('/lustre/aoc/students/jotter/dendro_catalogs/allbands_catalog.fits', format='fits', overwrite=True)

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
idx, d2d, d3d = dendro_coord.match_to_catalog_sky(Fb_coord)
#only qualify as a match if they are within 1"
matches = np.where(d2d.degree < 0.5*(1/3600))

for all_ind in matches[0]:
	Fb_data[idx[all_ind]]['D_ID'] = dendro_table[all_ind]['D_ID']

Fb_dend_table = join(dendro_table, Fb_data, keys='D_ID', join_type='left')


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
HC_data.remove_columns(['RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs'])

#doing the same for the RRS catalog
RRS_data = fits.getdata('/lustre/aoc/students/jotter/RRS2008_data.fit')

ra_RRS = RRS_data['_RAJ2000']*u.degree
dec_RRS = RRS_data['_DEJ2000']*u.degree
ID_RRS = RRS_data['Name'].strip()
ID_empty_RRS = np.full(len(ID_RRS),999999)

RRS_table = Table((ID_empty_RRS, ID_RRS, ra_RRS, dec_RRS), names=('HC_ID', 'RRS_ID', 'RA_RRS', 'DEC_RRS'))

#first need to match HC and RRS tables
HC_coord = SkyCoord(ra=HC_data['RA_HC'], dec=HC_data['DEC_HC'])
RRS_coord = SkyCoord(ra=RRS_table['RA_RRS'], dec=RRS_table['DEC_RRS'])
HC_idx, HC_d2d, HC_d3d = HC_coord.match_to_catalog_sky(RRS_coord)
#match within 0.5 arcsec
HC_match = np.where(HC_d2d.degree < 0.5*(1/3600))
#now loop through to add RRS IDs to HC table and vice versa
for match in HC_match[0]:
	HC_data[match]['RRS_ID'] = RRS_table[HC_idx[match]]['RRS_ID']
	RRS_table[HC_idx[match]]['HC_ID'] = HC_data[match]['HC_ID']

#ext_table is the combination of RA/decs/IDs from the HC and RRS catalogs
ext_table = join(HC_data, RRS_table,keys=('HC_ID','RRS_ID'), join_type='outer')


ext_table['RA_ext'] = np.zeros(len(ext_table))
ext_table['DEC_ext'] = np.zeros(len(ext_table))
ext_table['RA_ext'][np.where(ext_table['RA_RRS'].mask == False)] = ext_table['RA_RRS'][np.where(ext_table['RA_RRS'].mask == False)]
ext_table['RA_ext'][np.where(ext_table['RA_ext'] == 0)] = ext_table['RA_HC'][np.where(ext_table['RA_ext'] == 0)]
ext_table['DEC_ext'][np.where(ext_table['DEC_RRS'].mask == False)] = ext_table['DEC_RRS'][np.where(ext_table['DEC_RRS'].mask == False)]
ext_table['DEC_ext'][np.where(ext_table['DEC_ext'] == 0)] = ext_table['DEC_HC'][np.where(ext_table['DEC_ext'] == 0)]


ext_coord = SkyCoord(ra=ext_table['RA_ext']*u.degree, dec=ext_table['DEC_ext']*u.degree)
ext_table['D_ID'] = np.full(len(ext_coord),-1)

#matching the dendrogram catalog and the HC/RRS catalog (ext_table)
Fb_dend_coord = SkyCoord(ra=Fb_dend_table['RA']*u.degree, dec=Fb_dend_table['DEC']*u.degree)
idx, d2d, d3d = Fb_dend_coord.match_to_catalog_sky(ext_coord)
#only qualify as a match if they are within 1"
matched_inds = np.where(d2d.degree < (1/3600))

for all_ind in matched_inds[0]:
	ext_table[idx[all_ind]]['D_ID'] = Fb_dend_table[all_ind]['D_ID']

IR_matched_table = join(Fb_dend_table, ext_table, keys='D_ID', join_type='left', table_names=[names[ind -1], names[ind]])
IR_matched_table.write('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits',format='fits',overwrite=True)

#IR_source_table = join(dendro_table, ext_table, keys='D_ID', join_type='outer', table_names=[names[ind -1], names[ind]])
#ext_table.write('/lustre/aoc/students/jotter/dendro_catalogs/IR_only_catalog.fits',format='fits',overwrite=True)



'''
#SED text files
wavelengths = set()

SED_data = []
with open('/lustre/aoc/students/jotter/SED_list.csv', newline='') as myfile:  
    reader = csv.reader(myfile)
    i = 0
    for row in reader:    
        SED_data.append(row)
        for i in range(len(row)):
            if (i-2)%5 == 0 and i != 2:
                wavelengths.add(row[i-5]+' um flux')
                wavelengths.add(row[i-5]+' um error')
#table for SED values with columns corresponding to the object names and wavelength measurements
SED_tab = Table(names=(np.append(['HC ID'], np.append(['RRS ID'],np.array(sorted(list(wavelengths)))))))
#changing dtype of id rows
SED_tab['HC ID'] = SED_tab['HC ID'].astype(str)
SED_tab['RRS ID'] = SED_tab['RRS ID'].astype(str)
#loop through and fill in fluxes
for ind in range(len(SED_data)): 
    SED_tab.add_row()
    SED_tab['HC ID'][ind] = SED_data[ind][0]
    SED_tab['RRS ID'][ind] = SED_data[ind][1]
    for j in range(len(SED_data[ind])+1):
        if (j-2)%5 == 0 and j != 2:
            lam = SED_data[ind][j-5]
            flux = SED_data[ind][j-4]
            flux_err = SED_data[ind][j-3]
            SED_tab[lam+' um flux'][ind] = flux
            SED_tab[lam+' um error'][ind] = flux_err

full_matched_cat = join(matched_cat, SED_tab, join_type='left')

SED_tab.write('/lustre/aoc/students/jotter/SED_table_data.fits',format='fits')
full_matched_cat.write('/lustre/aoc/students/jotter/full_matched_cat.fits',format='fits')
'''    
