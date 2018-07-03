#This file matches the dendrogram catalogs with the HC2000 catalog and the RRS2008 catalog
from astropy.io import fits
from astropy import units as u
import numpy as np
from astropy.coordinates import SkyCoord, Angle
import matplotlib.pyplot as plt
from astropy.table import Table,hstack,vstack,join,Column
import csv
import glob

#first need to join together dendrogram catalogs
dendro_cats = sorted(glob.glob('/lustre/aoc/students/jotter/dendro_catalogs/*_dendro_catalog.fits'))
names = ['340GHz', '470GHz', 'B3', 'B6', 'B7']

dendro_table = Table(fits.getdata(dendro_cats[0]))
dendro_table['D_ID'] = np.arange(0,len(dendro_table))
dendro_table['RA'] = dendro_table['gauss_x '+names[0]]
dendro_table['DEC'] = dendro_table['gauss_y '+names[0]]

for ind in np.arange(1,len(dendro_cats)):
	print(ind)
	#need to match against all previous sources
	ra_all = dendro_table['RA']
	dec_all = dendro_table['DEC']
	#for j in np.arange(1,ind):
	#	ra_all = np.concatenate((ra_all, dendro_table['gauss_x '+names[j]]))
	#	dec_all = np.concatenate((dec_all, dendro_table['gauss_y '+names[j]]))
	coord_all = SkyCoord(ra=ra_all*u.degree, dec=dec_all*u.degree)
	#ABOVE IS THE PROBLEMATIC PART, NEED TO RETAIN POSITION IN TABLE
	#SOLUTION - create new columns for RA/DEC, just using the first one entered for that source

	cat = Table(fits.getdata(dendro_cats[ind]))
	cat['D_ID'] = np.full(len(cat),-1)
	ra2 = cat['gauss_x '+names[ind]]
	dec2 = cat['gauss_y '+names[ind]]
	coord2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
	idx, d2d, d3d = coord_all.match_to_catalog_sky(coord2)
	#idx is a list of indices for cat with the list index corresponding to the match in dendro_table
	match = np.where(d2d.degree < 0.1*(1/3600)) #matches are within 0.05 arcsec
	#match is indices of idx of matches
	#now loop through all sources
	for all_ind, cat_ind in enumerate(idx):
		if all_ind in match[0]:
			cat[idx[all_ind]]['D_ID'] = dendro_table[all_ind]['D_ID']
		
	dendro_table = join(dendro_table, cat, keys='D_ID', join_type='outer', table_names=[names[ind -1], names[ind]])
	#give unmatched sources a new, unique D_ID
	unmatched = np.where(dendro_table['D_ID'] == -1)
	for un in unmatched[0]:
		dendro_table[un]['D_ID'] = np.nanmax(np.concatenate((dendro_table['D_ID'], cat['D_ID'])))+1
		dendro_table[un]['RA'] = cat[un]['gauss_x ' + names[ind]]
		dendro_table[un]['DEC'] = cat[un]['gauss_y ' + names[ind]]
	
	
'''
HC_idx, HC_d2d, HC_d3d = HC_coord.match_to_catalog_sky(RRS_coord)
#match within 0.5 arcsec
HC_match = np.where(HC_d2d.degree < 0.5*(1/3600))
#now loop through to add RRS IDs to HC table and vice versa
for match in HC_match[0]:
	HC_table[match]['RRS ID'] = RRS_table[HC_idx[match]]['RRS ID']
	RRS_table[HC_idx[match]]['HC ID'] = HC_table[match]['HC ID']

#ext_table is the combination of RA/decs/IDs from the HC and RRS catalogs
ext_table = join(HC_table, RRS_table,keys=('HC ID', 'RRS ID'), join_type='outer')
'''
'''

my_file = '/lustre/aoc/students/jotter/catalog_pbcorr.fits'
my_cat = fits.getdata(my_file)

#ra/dec are given in degrees in the dendrogram catalog
ra_my = my_cat['x_cen']*u.degree
dec_my = my_cat['y_cen']*u.degree

my_coord = SkyCoord(ra=ra_my, dec=dec_my)

#loading in the HC catalog
col_locs = np.array([0,4,6,9,15,18,21,26,34,43,50,57])
widths = col_locs[1:] - col_locs[:-1]
HC_data = np.genfromtxt('/lustre/aoc/students/jotter/mr.table1.dat',skip_header=37,delimiter=widths,autostrip=True)
ID_HC = HC_data[:,0]
ID_HC = ID_HC.astype('int')
ID_empty_HC = np.full(len(HC_data[:,0]),'none',dtype='U10')
#converting from hours,minutes,seconds to degrees
ra_HC =(HC_data[:,1].astype(float)+HC_data[:,2].astype(float)/60+HC_data[:,3].astype(float)/3600)*(360/24)
ra_HC = ra_HC*u.degree
#degrees + arcmins and arcsecs
dec_HC = HC_data[:,4].astype(float)+(HC_data[:,5].astype(float)/60+HC_data[:,6].astype(float)/3600)*np.sign(HC_data[:,4].astype(float))
dec_HC = dec_HC*u.degree

HC_table = Table((ID_HC, ID_empty_HC, ra_HC, dec_HC, HC_data[:,7], HC_data[:,8], HC_data[:,9], HC_data[:,10]), names=('HC ID', 'RRS ID', 'RA_HC', 'DEC_HC', 'Kmag', 'Hmag', 'Kmag relative error', 'Hmag relative error'), dtype=('unicode_', 'unicode_', 'float', 'float', 'str', 'str', 'float', 'float'))

#doing the same for the RRS catalog
RRS_data = fits.getdata('/lustre/aoc/students/jotter/RRS2008_data.fit')

ra_RRS = RRS_data['_RAJ2000']*u.degree
dec_RRS = RRS_data['_DEJ2000']*u.degree
ID_RRS = RRS_data['Name'].strip()
ID_empty_RRS = np.full(len(ID_RRS),'none')

RRS_table = Table((ID_empty_RRS, ID_RRS, ra_RRS, dec_RRS), names=('HC ID', 'RRS ID', 'RA_RRS', 'DEC_RRS'))

#first need to match HC and RRS tables
HC_coord = SkyCoord(ra=HC_table['RA_HC'], dec=HC_table['DEC_HC'])
RRS_coord = SkyCoord(ra=RRS_table['RA_RRS'], dec=RRS_table['DEC_RRS'])
HC_idx, HC_d2d, HC_d3d = HC_coord.match_to_catalog_sky(RRS_coord)
#match within 0.5 arcsec
HC_match = np.where(HC_d2d.degree < 0.5*(1/3600))
#now loop through to add RRS IDs to HC table and vice versa
for match in HC_match[0]:
	HC_table[match]['RRS ID'] = RRS_table[HC_idx[match]]['RRS ID']
	RRS_table[HC_idx[match]]['HC ID'] = HC_table[match]['HC ID']

#ext_table is the combination of RA/decs/IDs from the HC and RRS catalogs
ext_table = join(HC_table, RRS_table,keys=('HC ID', 'RRS ID'), join_type='outer')
RA_ext = np.
ext_coord = SkyCoord(ra=ext_table['RA_ext'], dec=ext_table['DEC_ext'])


#matching the dendrogram catalog and the HC/RRS catalog (ext_table)
idx, d2d, d3d = my_coord.match_to_catalog_sky(ext_coord)
#idx is a list of indices for ext_coord with the list index corresponding to the match in my_coord
matched_tab = ext_coord[idx]

#only qualify as a match if they are within 1"
matched_inds = np.where(d2d.degree < (1/3600))
#indices of unmatched sources
unmatched_inds = np.delete(np.arange(len(d2d)),matched_inds)

#combining matches in dendrogram catalogue and HC/RRS coordinates
d2d_tab = Column(d2d.degree, name='match_dist')
matched_tab_coord = ext_coord[idx[matched_inds]]
matched_cat = hstack((Table(my_cat[matched_inds]), Table(ext_table[idx[matched_inds]])),join_type='outer')
matched_cat.add_column(d2d_tab[matched_inds], name='match_dist')
matched_cat['match_dist'].unit = 'degree'

#catalogue for unmatched sources
unmatched_cat = Table(my_cat[unmatched_inds])
bad_sources = [4,11,12,17,21,22,23,25,28,29,36,60,61,70,74,78,98,100]
mask = np.in1d(unmatched_cat['_idx'], bad_sources)
unmatched_cat.remove_rows(np.where(mask))

#unmatched_cat.write('/lustre/aoc/students/jotter/unmatched_cat.fits',format='fits')

#not worrying about the SED text files for now


plt.scatter(ext_table['RA_HC'], ext_table['DEC_HC'], c='g',alpha=0.3)
plt.scatter(ra_my, dec_my, c='r',alpha=0.3)
plt.scatter(ext_table['RA_RRS'], ext_table['DEC_RRS'], c='k', alpha=0.6)
plt.xlabel('RA (degrees)')
plt.ylabel('dec (degrees)')
plt.show()
'''

'''
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
