from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import gaussian_kde, ks_2samp


def plot_MST_hist(coord_list, labels, filename, colors=['tab:blue', 'tab:orange'], KS_test=False):

    fig = plt.figure(figsize=(10,6))
    dists = []

    means = []
    errs = []

    for i, tab_coord in enumerate(coord_list):
        graph = []

        for src_coord in tab_coord:
            graph.append(src_coord.separation(tab_coord).arcsecond)

        graph = csr_matrix(graph)
        mst_graph = minimum_spanning_tree(graph)
        mst_graph = mst_graph.toarray()

        sep_ind = np.where(mst_graph != 0)
        sep_arr = mst_graph[sep_ind]

        dists.append(sep_arr)

        #bins = np.arange(0,15.1,1.5)
        bins = np.arange(0,20,1.5)
        print(f'Mean separation: {np.mean(sep_arr)} arcseconds for {labels[i]}')
        print(f'Median separation: {np.median(sep_arr)} arcseconds for {labels[i]}')
        print(f'standard dev on sep: {np.std(sep_arr)} arcseconds for {labels[i]}')
        print(f'standard error on sep: {np.std(sep_arr)/np.sqrt(len(sep_arr))} arcseconds for {labels[i]}')

        means.append(np.mean(sep_arr))
        errs.append(np.std(sep_arr)/np.sqrt(len(sep_arr)))

        sep_au = np.mean(sep_arr) * 400
        sep_au_med = np.median(sep_arr) * 400
        sep_err_au = np.std(sep_arr)/np.sqrt(len(sep_arr)) * 400
        sep_3d = sep_au * 1.3
        sep_3d_med = sep_au_med * 1.3
        sep_3d_err = sep_err_au * 1.3
        
        src_density = 3 / (4*np.pi*((sep_3d*u.AU).to(u.pc))**3)
        src_density_med = 3 / (4*np.pi*((sep_3d_med*u.AU).to(u.pc))**3)
        src_density_err = 9 / (4*np.pi*((sep_3d*u.AU).to(u.pc))**4) * (sep_3d_err*u.AU).to(u.pc)

        print(f'Source density: {src_density} pm {src_density_err} pc^-3')
        print(f'Source density (median): {src_density_med} pm {src_density_err} pc^-3')
        
        
        sep_grid = np.linspace(0,20,100)
        sep_kde = gaussian_kde(sep_arr, bw_method='scott')
        print(sep_kde.factor)
        sep_pdf = sep_kde.evaluate(sep_grid)                                                                                                                                                         #to normalize, multiply by bin width and number of points
        sep_est = sep_pdf*len(sep_arr)*(bins[1]-bins[0])

        plt.hist(sep_arr, bins=bins, label=labels[i], alpha=0.4, color=colors[i], edgecolor='black')
        plt.plot(sep_grid, sep_est, label=labels[i]+' KDE', color=colors[i])

    if KS_test == True:
        Dval, pval = ks_2samp(dists[0], dists[1])
        print(f'{pval} p-value of KS test between first two coord lists')

    discrep = np.abs(means[1] - means[0])
    errtot = np.sqrt(errs[0]**2 + errs[1]**2)

    print(f't value: {discrep/errtot}')
    
    plt.legend()
    plt.xlabel('Separation (arcsecond)')
    plt.ylabel('Number')
    plt.savefig(filename, bbox_inches='tight')




tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
tab2 = Table.read('/home/jotter/nrao/tables/eis_coord_table.fits')

IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
nonIR_src = np.setdiff1d(tab['ID'], IR_tab['ID'])
nonIR_ind = [np.where(tab['ID']==d_id)[0][0] for d_id in nonIR_src]

#nonIR_ind = nonIR_ind[:-1] #to remove source 124 from OMC1 sample (outlier)
#print(nonIR_ind)

IR_src = IR_tab['ID']
IR_ind = [np.where(tab['ID']==d_id)[0][0] for d_id in IR_src]

tab_coord = SkyCoord(ra=tab['RA_B3'], dec=tab['DEC_B3'], unit=u.degree)
tab_coord_IR = SkyCoord(ra=tab['RA_B3'][IR_ind], dec=tab['DEC_B3'][IR_ind], unit=u.degree)
tab_coord_nonIR = SkyCoord(ra=tab['RA_B3'][nonIR_ind], dec=tab['DEC_B3'][nonIR_ind], unit=u.degree)


tab_coord2 = SkyCoord(ra=tab2['RA'], dec=tab2['DEC'], unit=u.degree)

coord_list = [tab_coord_IR, tab_coord2]

plot_MST_hist(coord_list, labels=['Band 3 ONC', 'E18'], filename='/home/jotter/nrao/plots/separation_hist_ONC.pdf', KS_test=True, colors=['tab:blue', 'tab:green'])
#plot_MST_hist(coord_list, labels=['Band 3 OMC1', 'Band 3 ONC', 'E18'], filename='/home/jotter/nrao/plots/separation_hist_all.pdf', colors=['tab:green', 'tab:blue', 'tab:orange'])

coord_list = [tab_coord_nonIR, tab_coord2]

plot_MST_hist(coord_list, labels=['Band 3 OMC1', 'E18'], filename='/home/jotter/nrao/plots/separation_hist_OMC1.pdf', KS_test=True, colors=['tab:orange','tab:green'])

#coord_list = [tab_coord_nonIR, tab_coord_IR]

#plot_MST_hist(coord_list, labels=['Band 3 OMC1', 'Band 3 ONC'], filename='/home/jotter/nrao/plots/separation_hist_oncomc1.pdf')
