from astropy.table import Table
from astropy.io import ascii
import numpy as np
from astropy.coordinates import SkyCoord
import regions
import astropy.units as u

def RA_to_deg(HH, MM, SS):
    return (HH + MM/60 + SS/(60*60))*(360/24)

def DEC_to_deg(DD, MM, SS):
    return DD + (MM/60 + SS/(60*60))*np.sign(DD)

eisner = ascii.read('../tables/eisner_tbl.txt', format='tab')
reg_file = '/users/jotter/summer_research_2018/final_regs/eisner_regs.reg'

with open(reg_file, 'w') as fh:
    fh.write("fk5\n")
    rad = (0.1*u.arcsec).to(u.deg).value
    for row in range(len(eisner)):
        ra = (RA_to_deg(float(eisner[row]['alpha'][0:2]), float(eisner[row]['alpha'][2:4]), float(eisner[row]['alpha'][4:])))
        dec = (DEC_to_deg(float(eisner[row]['delta'][0:3]), float(eisner[row]['delta'][3:5]), float(eisner[row]['delta'][5:])))
        #reg = regions.CircleSkyRegion(center=SkyCoord(ra*u.degree, dec*u.degree), radius=rad, meta={'text':str(eisner[row]['ID'])})
        fh.write('circle({x_cen}, {y_cen}, {radius}) #text={{{ID}}}\n'.format(x_cen=ra, y_cen=dec, radius=rad, ID=str(eisner[row]['ID'])+' eis'))

        
data = Table.read('../tables/r0.5_catalog_conv_bgfitted_add_final3.fits')
reg_file2 = '/users/jotter/summer_research_2018/final_regs/final_cat_regs.reg'

with open(reg_file2, 'w') as fh:
    fh.write("fk5\n")
    rad = (0.1*u.arcsec).to(u.deg).value
    for row in range(len(data)):
        ra = data['RA_B3'][row]
        dec = data['DEC_B3'][row]
        #reg = regions.CircleSkyRegion(center=SkyCoord(ra*u.degree, dec*u.degree), radius=rad, meta={'text':str(eisner[row]['ID'])})
        fh.write('circle({x_cen}, {y_cen}, {radius}) #text={{{ID}}}\n'.format(x_cen=ra, y_cen=dec, radius=rad, ID=str(data['D_ID'][row])+' r0.5'))


reg_file3 = '/users/jotter/summer_research_2018/final_regs/final_cat_regs_B6.reg'

data_B6 = data[np.where(np.isnan(data['ap_flux_B6']) == False)[0]]
with open(reg_file3, 'w') as fh:
    fh.write("fk5\n")
    rad = (0.1*u.arcsec).to(u.deg).value
    for row in range(len(data_B6)):
        ra = data_B6['RA_B3'][row]
        dec = data_B6['DEC_B3'][row]
        fh.write('circle({x_cen}, {y_cen}, {radius}) #text={{{ID}}}\n'.format(x_cen=ra, y_cen=dec, radius=rad, ID=str(data['D_ID'][row])+' r0.5'))

reg_file4 = '/users/jotter/summer_research_2018/final_regs/final_cat_regs_B7.reg'

data_B7 = data[np.where(np.isnan(data['ap_flux_B7']) == False)[0]]
with open(reg_file4, 'w') as fh:
    fh.write("fk5\n")
    rad = (0.1*u.arcsec).to(u.deg).value
    for row in range(len(data_B7)):
        ra = data_B7['RA_B3'][row]
        dec = data_B7['DEC_B3'][row]
        fh.write('circle({x_cen}, {y_cen}, {radius}) #text={{{ID}}}\n'.format(x_cen=ra, y_cen=dec, radius=rad, ID=str(data['D_ID'][row])+' r0.5'))
