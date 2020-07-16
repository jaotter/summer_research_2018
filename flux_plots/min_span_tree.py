from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import gaussian_kde


def plot_MST_hist(coord_list, labels, filename, colors=['tab:blue', 'tab:orange']):

    fig = plt.figure(figsize=(10,6))
    for i, tab_coord in enumerate(coord_list):
        
        graph = []

        for src_coord in tab_coord:
            graph.append(src_coord.separation(tab_coord).arcsecond)

        graph = csr_matrix(graph)
        mst_graph = minimum_spanning_tree(graph)
        mst_graph = mst_graph.toarray()

        sep_ind = np.where(mst_graph != 0)
        sep_arr = mst_graph[sep_ind]

        bins = np.arange(0,15.1,1.5)
        print(f'Mean separation: {np.mean(sep_arr)} arcseconds for {labels[i]}')

        sep_grid = np.linspace(0,14,100)
        sep_kde = gaussian_kde(sep_arr)
        sep_pdf = sep_kde.evaluate(sep_grid)                                                                                                                                                         #to normalize, multiply by bin width and number of points
        sep_est = sep_pdf*len(sep_arr)*(bins[1]-bins[0])

        plt.hist(sep_arr, bins=bins, label=labels[i], alpha=0.5, color=colors[i], edgecolor='black')
        plt.plot(sep_grid, sep_est, label=labels[i]+' KDE', color=colors[i])

    plt.legend()
    plt.xlabel('Separation (arcsecond)')
    plt.ylabel('Number')
    plt.savefig(filename)




tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_apr20.fits')

tab2 = Table.read('/home/jotter/nrao/tables/eis_coord_table.fits')

tab_coord = SkyCoord(ra=tab['RA_B3'], dec=tab['DEC_B3'], unit=u.degree)
tab_coord2 = SkyCoord(ra=tab2['RA'], dec=tab2['DEC'], unit=u.degree)

coord_list = [tab_coord, tab_coord2]

plot_MST_hist(coord_list, labels=('Band 3', 'E18'), filename='/home/jotter/nrao/plots/separation_hist.pdf')
