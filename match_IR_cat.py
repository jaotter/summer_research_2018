from astropy.table import Table
from astropy.coordinates import SkyCoord

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

MLLA = Table.read('/home/jotter/nrao/tables/MLLA_02_IR.fit')
A11 = Table.read('/home/jotter/nrao/tables/A11_IR.fit')
data = Table.read('/home/jotter/nrao/tables/r0.5_catalog_conv_bgfitted_add_final3_ann2.fits')

MLLA_coord = SkyCoord(ra=MLLA['RAJ2000'], dec=MLLA['DEJ2000'], unit=u.degree)
A11_coord = SkyCoord(ra=A11['RAJ2000'], dec=A11['DEJ2000'], unit=u.degree)
B3_coord = SkyCoord(ra=data['RA_B3'], dec=data['DEC_B3'], unit=u.degree)

idx_MLLA, d2d_MLLA, d3d = B3_coord.match_to_catalog_sky(MLLA_coord)
not_matches_MLLA = np.where(d2d_MLLA > 0.5*u.arcsecond)

MLLA_nm = MLLA[idx_MLLA[not_matches_MLLA[0]]]


idx_A11, d2d_A11, d3d = B3_coord.match_to_catalog_sky(A11_coord)
not_matches_A11 = np.where(d2d_A11 > 0.5*u.arcsecond)

A11_nm = A11[idx_A11[not_matches_A11[0]]]


#compare non matches to our FOV
