from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord

#import regions
import numpy as np

b3_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_feb21_ulim.fits')
MLLA = Table.read('/home/jotter/nrao/tables/MLLA_02_IR.fit')

b3_center = SkyCoord(ra=83.8104625, dec=-5.37515556, unit=u.degree)
b3_radius = 137.6*u.arcsecond/2

#b3_region = CircleSkyRegion(b3_center, b3_radius)

mlla_coord = SkyCoord(ra=MLLA['RAJ2000'], dec=MLLA['DEJ2000'], unit=u.degree)

seps = mlla_coord.separation(b3_center)
within_fov_ind = np.where(seps < b3_radius)

print(len(within_fov_ind[0]))
print(len(within_fov_ind[0]) - 67)


MLLA_match = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_feb21_full.fits')
mlla_match_coord = SkyCoord(ra=MLLA_match['RAJ2000'], dec=MLLA_match['DEJ2000'], unit=u.degree)

seps = mlla_match_coord.separation(b3_center)
within_fov_ind = np.where(seps < b3_radius)

print(len(within_fov_ind[0]))
print(len(mlla_match_coord))
