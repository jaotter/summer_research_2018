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
    

#data = Table.read('/home/jotter/nrao/tables/HC2000.fit')
#reg_file = '/home/jotter/nrao/images/regions/HC2000_IR.reg'
#print(data.info)
#data = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full.fits')
#reg_file = '/home/jotter/nrao/images/regions/IR_matches_may21.reg'

#data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_new_det_mar21.fits')
#reg_file = '/home/jotter/nrao/images/regions/r0.5_new_det_mar21.reg'

#data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim.fits')
#reg_file = '/home/jotter/nrao/images/regions/r0.5_poorfit_may21.reg'
#ind = np.array([113,114,124,125])

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/OMC1_r0.5_may21.fits')
reg_file = '/home/jotter/nrao/images/regions/OMC1_may21.reg'


#data = Table.read('/home/jotter/nrao/tables/MLLA_02_IR.fit')
#reg_file = '/home/jotter/nrao/images/regions/MLLA_full.reg'

#data = Table.read('/home/jotter/nrao/tables/COUP_srclist.fits')
#reg_file = '/home/jotter/nrao/images/regions/COUP_full.reg'

#data = Table.read('/home/jotter/nrao/tables/robberto2010.fits')
#reg_file = '/home/jotter/nrao/images/regions/robberto2010.reg'

#coord = SkyCoord(ra=data['RAh'], dec=data['DEdeg'], unit=(u.hourangle,u.degree))

#write_reg_file(coord.ra.deg, coord.dec.deg, reg_file, name='R10', IDs=data['ID'])

#data = Table.read('/home/jotter/nrao/tables/eis_nondet_full.fits')
#reg_file = '/home/jotter/nrao/images/regions/E18_nondet.reg'


#data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_mar21_ulim.fits')
#lowsnr = np.where(data['SNR_B3'] < 5)
#nob6 = np.where(np.isnan(data['SNR_B6']) == True)
#nob7 = np.where(np.isnan(data['SNR_B6']) == True)
#nob6b7 = np.intersect1d(nob6, nob7)
#b3_lowsnr = np.intersect1d(lowsnr, nob6b7)
#reg_file = '/home/jotter/nrao/images/regions/b3_lowsnr.reg'

#data = Table.read('/home/jotter/nrao/tables/r0.5_b3_catalog_full_may21.fits')
#reg_file = '/home/jotter/nrao/images/regions/b3_full_b3seq_may21.reg'

#data = Table.read('/home/jotter/nrao/summer_research_2018/tables/COUP_may21_nondet_nonIR.fits')
#reg_file = '/home/jotter/nrao/images/regions/coup_nondet_nonIR.reg'

data = Table.read('/home/jotter/nrao/summer_research_2018/tables/Forbrich16_may21_nondet_nonIR.fits')
reg_file = '/home/jotter/nrao/images/regions/fb16_nondet_nonIR.reg'


write_reg_file(data['RAJ2000'], data['DEJ2000'], reg_file, name='', IDs=data['Seq'])
