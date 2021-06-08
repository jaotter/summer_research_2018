from astropy.table import Table, Column
from astropy.io import fits, ascii

import numpy as np

def RA_to_deg(HH, MM, SS):
    return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
    return DD + (MM/60 + SS/(60*60))*np.sign(DD)

def reg_to_table(reg_file, name=None):
    reg = open(reg_file, mode='r')
    ra_arr = []
    dec_arr = []
    reg_lines = reg.readlines()
    for line_ind in np.arange(3,len(reg_lines)):
        line = reg_lines[line_ind]
        split_par = line.split('(')[1]
        split_comma = split_par.split(',')
        ra = split_comma[0]
        dec = split_comma[1]
        rasplit = ra.split(':')
        decsplit = dec.split(':')

        ra_arr.append(RA_to_deg(float(rasplit[0]), float(rasplit[1]), float(rasplit[2])))
        dec_arr.append(DEC_to_deg(float(decsplit[0]), float(decsplit[1]), float(decsplit[2])))
    reg.close()

    tab = Table([np.arange(len(ra_arr)), ra_arr, dec_arr], names=(f'Seq{"_"+name if name else ""}', f'RA{"_"+name if name else ""}', f'DEC{"_"+name if name else ""}'))
    return tab

#B3_reg = '/home/jotter/nrao/images/dendro_regions/B3_cleaned_minval6e-05_mindel9e-05_npix10.reg'
#B3_tab = reg_to_table(B3_reg, name='B3')
#B3_tab.write('/home/jotter/nrao/tables/B3_dendro_table_cleaned.fits', format='fits')

#B6_reg = '/home/jotter/nrao/images/dendro_regions/B6_cleaned_minval0.00018_mindel0.00035_npix10.reg'
#B6_tab = reg_to_table(B6_reg, name='B6')
#B6_tab.write('/home/jotter/nrao/tables/B6_dendro_table_cleaned.fits', format='fits')

#B7_reg = '/home/jotter/nrao/images/dendro_regions/B7_cleaned_minval0.00065_mindel0.0012_npix10.reg'
#B7_tab = reg_to_table(B7_reg, name='B7')
#B7_tab.write('/home/jotter/nrao/tables/B7_dendro_table_cleaned.fits', format='fits')

B3_reg = '/home/jotter/nrao/images/regions/b3_huge_visual_ID.reg'
B3_tab = reg_to_table(B3_reg, name='B3')
B3_tab.write('/home/jotter/nrao/tables/B3_huge_visual_ID.fits', format='fits')
