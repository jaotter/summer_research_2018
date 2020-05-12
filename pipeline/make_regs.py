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

def write_reg_file(ra, dec, filename, name='', IDs=None):
    with open(filename, 'w') as fh:
        fh.write("fk5\n")
        rad = (0.1*u.arcsec).to(u.deg).value
        for ind in range(len(ra)):
            ID_str = f'{name} {IDs[ind] if IDs is not None else ""}'
            fh.write('circle({x_cen}, {y_cen}, {radius}) #text={{{ID}}}\n'.format(x_cen=ra[ind], y_cen=dec[ind], radius=rad, ID=ID_str))
    

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_apr20.fits')
reg_file = '/home/jotter/nrao/images/regions/r0.5_catalog_apr20.reg'

write_reg_file(data['RA_B3'], data['DEC_B3'], reg_file, name='r0.5', IDs=data['D_ID'])
