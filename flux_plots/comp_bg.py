from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from flux_plots import get_ind

data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits')
ind = get_ind(['B3', 'B6', 'B7_hr', 'B7_lr'])



plt.scatter(data['D_ID'][ind], data['bg_mean_B7_lr'][ind], marker='o', label='lr')
plt.scatter(data['D_ID'][ind], data['bg_mean_B7_hr'][ind], marker='o', label='hr')
plt.xlabel('source number')
plt.ylabel('LR background/HR background')
plt.show()
