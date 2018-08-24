import numpy as np
from astropy.table import Table
from functools import reduce
from flux_plots import multi_img_SED
import matplotlib.pyplot as plt

data = Table.read('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits')
RA_names = ['gauss_x_B3', 'gauss_x_B6', 'gauss_x_B7_hr']
inds = []
for fn in RA_names:
	ind = np.where(np.isnan(data[fn]) == False)
	inds.append(ind)

detected = reduce(np.intersect1d, inds)

for src in detected:
	print(src)
	srcID = data['D_ID'][src]
	multi_img_SED(srcID, 'r-2.clean0.1mJy.500klplus.deepmask', 'r-2.clean0.1mJy.500klplus.deepmask', 'r-2.clean0.1mJy.500klplus.deepmask', 'r-2.clean0.1mJy.500klplus.deepmask', flux_type='ap')

