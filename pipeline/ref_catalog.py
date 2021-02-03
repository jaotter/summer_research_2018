from astropy.table import Table, Column
from astropy.io import fits, ascii

import numpy as np

def RA_to_deg(HH, MM, SS):
    return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
    return DD + (MM/60 + SS/(60*60))*np.sign(DD)
        

table = Table.read('/home/jotter/nrao/summer_research_2018/tables/B3_visual_ID_matched.fits')

ref_table = Table((table['Seq_B3'], table['RA_B3_1'], table['DEC_B3_1'], table['D_ID']))
ref_table.add_column(Column(np.logical_not(np.isnan(table['ap_flux_B3'])), name='B3_detect'))
ref_table.add_column(Column(np.logical_not(np.isnan(table['ap_flux_B6'])), name='B6_detect'))
ref_table.add_column(Column(np.logical_not(np.isnan(table['ap_flux_B7'])), name='B7_detect'))
ref_table.rename_columns(['RA_B3_1', 'DEC_B3_1'], ['RA_B3', 'DEC_B3'])
ref_table['B3_detect'][ref_table['D_ID'] > 1000] = True

#ind1 = np.where(ref_table['D_ID'] == 26)[0]
#ref_table['B6_detect'][ind1] = True

#ind2 = np.where(ref_table['D_ID'] == 52)[0]
#ref_table['B6_detect'][ind2] = True

#for i  in range(len(ras)):
#    new_id = ref_table['D_ID'][-1]+1
#    new_ra = ras[i]
#    new_dec = decs[i]
#    ref_table.add_row([new_id, new_ra, new_dec, True, False, False])

ref_table.write('../tables/ref_catalog_feb21.fits', overwrite=True)
