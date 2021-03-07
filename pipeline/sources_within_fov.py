from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord

import regions
import numpy as np

b3_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_feb21_ulim.fits')
MLLA = Table.read('/home/jotter/nrao/tables/MLLA_02_IR.fit')

b3_center = SkyCoord(ra=83.8104625, dec=-5.37515556, unit=u.degree)
