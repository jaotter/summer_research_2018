from astropy.table import Table
import numpy as np

tab = Table.read('/home/jotter/nrao/tables/array_locs_B7.txt', format='ascii')

north = tab['North']
east = tab['East']
elevation = tab['Elevation']


max_dist_total = 0
min_dist_total = 100000

max_antennae = []
min_antennae = []

for i in range(len(north)):
    e_dists = east - east[i] 
    n_dists = north - north[i]
    ele_dists = elevation - elevation[i]
    total_dists = np.sqrt(e_dists**2 + n_dists**2 + ele_dists**2)
    min_dist = np.min(total_dists[total_dists > 0])
    max_dist = np.max(total_dists)

    if max_dist > max_dist_total:
        max_dist_total = max_dist
        max_antennae = [i, np.where(total_dists == max_dist)[0][0]]
    if min_dist < min_dist_total:
        min_dist_total = min_dist
        min_antennae = [i, np.where(total_dists == min_dist)[0][0]]


print(f'Maximum baseline distance: {max_dist_total} m for antennae {max_antennae}')
print(f'Minimum baseline distance: {min_dist_total} m for antennae {min_antennae}')
