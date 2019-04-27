from astropy.table import Table, Column
from astropy.io import fits, ascii

import numpy as np

def RA_to_deg(HH, MM, SS):
    return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
    return DD + (MM/60 + SS/(60*60))*np.sign(DD)
        

table = Table.read('../tables/r0.5_catalog_conv_bgfitted_apflux_final.fits')
eis_table = Table.read('../tables/eis_coord_table.fits')

ref_table = Table((table['D_ID'], table['RA_B3'], table['DEC_B3']))
ref_table.add_column(Column(np.logical_not(np.isnan(table['ap_flux_B3'])), name='B3_detect'))
ref_table.add_column(Column(np.logical_not(np.isnan(table['ap_flux_B6'])), name='B6_detect'))
ref_table.add_column(Column(np.logical_not(np.isnan(table['ap_flux_B7'])), name='B7_detect'))

ind1 = np.where(ref_table['D_ID'] == 26)[0]
ref_table['B6_detect'][ind1] = True

ind2 = np.where(ref_table['D_ID'] == 52)[0]
ref_table['B6_detect'][ind2] = True

ras = [RA_to_deg(5,35,14.0589), RA_to_deg(5,35,16.3130), RA_to_deg(5,35,15.889), RA_to_deg(5,35,16.1474), RA_to_deg(5,35,16.0796), RA_to_deg(5,35,15.6387)]
decs = [DEC_to_deg(-5,22,5.661), DEC_to_deg(-5,22,21.525), DEC_to_deg(-5,22,33.19), DEC_to_deg(-5,22,55.313), DEC_to_deg(-5,22,54.345), DEC_to_deg(-5,22,56.453)]
#corresponding eisner IDs: none, 163-222, HC447, HC393, HC401, HC389

for i  in range(len(ras)):
    new_id = ref_table['D_ID'][-1]+1
    new_ra = ras[i]
    new_dec = decs[i]
    ref_table.add_row([new_id, new_ra, new_dec, True, False, False])

ref_table.write('../tables/ref_catalog_added.fits', overwrite=True)
