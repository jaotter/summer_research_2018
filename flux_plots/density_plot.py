from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.spatial import KDTree
from matplotlib.patches import Circle
import astropy.wcs.utils as wcs_utils
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def make_density_map(npix, n_neighbor, savepath=None, ax_input=None, pos=111,
                     sources='all', vbounds=(0,100), plot_pb=False):
    if sources == 'all':
        tab = Table.read('tables/r0.5_catalog_bgfit_mar21_ulim.fits')
    if sources == 'IR':
        tab = Table.read('tables/IR_matches_MLLA_nov20_full.fits')
    if sources == 'nonIR' or sources == 'all_contour':
        full_tab = Table.read('tables/r0.5_catalog_bgfit_mar21_ulim.fits')
        tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
    if sources == 'IR':
        tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
    if sources == 'nonIR' or sources == 'all_contour':
        full_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')
        IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
        IR_tab = Table.read('tables/IR_matches_MLLA_nov20_full.fits')
        IR_tab = Table.read('tables/IR_matches_MLLA_mar21_full.fits')

        nonIR_src = np.setdiff1d(full_tab['ID'], IR_tab['ID'])
        nonIR_ind = [np.where(full_tab['ID']==d_id)[0][0] for d_id in nonIR_src]
        tab = full_tab[nonIR_ind]

    if sources == 'all_contour':
        tab_omc1 = tab
        #tab = Table.read('tables/IR_matches_MLLA_nov20_full.fits')
        tab = Table.read('tables/IR_matches_MLLA_mar21_full.fits')
    b3_fl = fits.open('images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits')
    b3_wcs = WCS(b3_fl[0].header).celestial
    b3_data = b3_fl[0].data.squeeze()

    b3_xmax = len(b3_data)
    b3_ymax = len(b3_data[0])

    orig_npix = b3_fl[0].header['NAXIS1']

    #creating new wcs for density map
    zero_coord = b3_wcs.all_pix2world(0,0,0)

    new_header = b3_fl[0].header
    new_header['CDELT1'] = b3_wcs.wcs.cdelt[0]*(b3_xmax/npix)
    new_header['CDELT2'] = b3_wcs.wcs.cdelt[1]*(b3_ymax/npix)
    new_header['CRPIX1'] = 1
    new_header['CRPIX2'] = 1
    new_header['CRVAL1'] = float(zero_coord[0])
    new_header['CRVAL2'] = float(zero_coord[1])
    new_header['NAXIS1'] = npix
    new_header['NAXIS2'] = npix

    new_wcs = WCS(new_header).celestial

    #density_map = np.zeros((npix,npix))
    xx, yy = np.mgrid[0:npix, 0:npix] #grid with density map pixel coords

    pix_points = np.array(list(zip(xx.ravel()+0.5,yy.ravel()+0.5)))

    src_coords = SkyCoord(tab['RA_B3'], tab['DEC_B3'], unit=u.degree)
    src_pix_coords = np.array(new_wcs.all_world2pix(src_coords.ra, src_coords.dec, 0))
    src_pix_coords = np.array([src_pix_coords[1], src_pix_coords[0]])
    src_map_coords = src_pix_coords
    src_pix_coords = np.transpose(src_pix_coords)

    tree = KDTree(src_pix_coords) #KDTree has source positions in map pix coords

    full_dists, full_inds = tree.query(pix_points, k=n_neighbor) #each entry of full_dists has nearest neighbors up to the nth nearest

    nth_dists_pix = [dist[n_neighbor-1] for dist in full_dists]

    #convert to arcseconds
    pixel_scales = wcs_utils.proj_plane_pixel_scales(b3_wcs)
    if pixel_scales[0] != pixel_scales[1]:
        print('WARNING: X/Y PIXEL SCALES DIFFERENT, DISTANCES WILL BE WRONG')

    pixel_scale = (pixel_scales[0]*u.degree).to(u.arcsecond)
    pixel_scale = (orig_npix / npix) * pixel_scale
    nth_dists = (nth_dists_pix * pixel_scale).value

    density_map = np.reshape(nth_dists, (npix,npix))


    if sources == 'all_contour':
        src_coords = SkyCoord(tab_omc1['RA_B3'], tab_omc1['DEC_B3'], unit=u.degree)
        src_pix_coords = np.array(new_wcs.all_world2pix(src_coords.ra, src_coords.dec, 0))
        src_pix_coords = np.array([src_pix_coords[1], src_pix_coords[0]])
        src_map_coords = src_pix_coords
        src_pix_coords = np.transpose(src_pix_coords)

        tree = KDTree(src_pix_coords) #KDTree has source positions in map pix coords

        full_dists, full_inds = tree.query(pix_points, k=n_neighbor) #each entry of full_dists has nearest neighbors up to the nth nearest

        nth_dists_pix = [dist[n_neighbor-1] for dist in full_dists]

        #convert to arcseconds
        nth_dists = (nth_dists_pix * pixel_scale).value

        density_map_omc1 = np.reshape(nth_dists, (npix,npix))

    #make plot
    if ax_input is None:
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=new_wcs)
    else:
        ax = ax_input
    im = ax.imshow(density_map, origin='lower', cmap=plt.cm.inferno_r, vmin=vbounds[0], vmax=vbounds[1])
    #plt.colorbar(im)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if plot_pb == True:
        B3_coord = SkyCoord(ra=83.8104625, dec=-5.37515556, unit=u.degree)
        B3_coord_pix = new_wcs.all_world2pix(B3_coord.ra, B3_coord.dec, 0)
        B3_radius = 60*u.arcsecond/2
        B3_coord_R = SkyCoord(ra=B3_coord.ra, dec=B3_coord.dec+B3_radius)
        B3_R_pix = new_wcs.all_world2pix(B3_coord_R.ra, B3_coord_R.dec,0)
        B3_radius_pix = B3_coord_pix[1] - B3_R_pix[1]
        B3_circ = Circle((B3_coord_pix[0], B3_coord_pix[1]), radius=B3_radius_pix, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='--')
        ax.add_patch(B3_circ)

        B3_fov = 137.6*u.arcsecond/2
        B3_coord_R = SkyCoord(ra=B3_coord.ra, dec=B3_coord.dec+B3_fov)
        B3_R_pix = new_wcs.all_world2pix(B3_coord_R.ra, B3_coord_R.dec,0)
        B3_radius_pix = B3_coord_pix[1] - B3_R_pix[1]
        B3_circ = Circle((B3_coord_pix[0], B3_coord_pix[1]), radius=B3_radius_pix, transform=ax.transData, color='blue', linewidth=1, fill=False, linestyle='-')
        ax.add_patch(B3_circ)

    if sources != 'all_contour':
        ax.scatter(src_map_coords[1], src_map_coords[0], marker='*')
    elif sources == 'all_contour':
        #c_levels = np.linspace(vbounds[0], vbounds[1]/2, 6)
        c_levels = [0.02,0.04,0.06,0.08]
        CS = ax.contour(density_map_omc1, c_levels, colors='tab:green')
        ax.clabel(CS, inline=True, fontsize=12, colors='tab:green')
        #plt.colorbar(im)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')

    return ax, im

npix=64
nth_neighbor=5
savepath = f'density_map_npix{npix}_{nth_neighbor}neighbor_3panel.pdf'

b3_fl = fits.open('images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
b3_wcs = WCS(b3_fl[0].header).celestial

zero_coord = b3_wcs.all_pix2world(0,0,0)
b3_xmax=len(b3_fl[0].data)
b3_ymax=b3_xmax
npix=64
new_header = b3_fl[0].header
new_header['CDELT1'] = b3_wcs.wcs.cdelt[0]*(b3_xmax/npix)
new_header['CDELT2'] = b3_wcs.wcs.cdelt[1]*(b3_ymax/npix)
new_header['CRPIX1'] = 1
new_header['CRPIX2'] = 1
new_header['CRVAL1'] = float(zero_coord[0])
new_header['CRVAL2'] = float(zero_coord[1])
new_header['NAXIS1'] = npix
new_header['NAXIS2'] = npix

new_wcs = WCS(new_header).celestial


#fig, [IR_ax, nonIR_ax, all_ax] = plt.subplots(1,3, subplot_kw={'projection':b3_wcs}, figsize=(15,4), constrained_layout=True)
fig = plt.figure(figsize=(15.75,5), constrained_layout=False)
gs = gridspec.GridSpec(1,3, figure=fig, wspace=0.05)
IR_ax = fig.add_subplot(gs[0], projection=new_wcs)
nonIR_ax = fig.add_subplot(gs[1], projection=new_wcs)
all_ax = fig.add_subplot(gs[2], projection=new_wcs)

IR_ax, im = make_density_map(npix,nth_neighbor,savepath=None,sources='IR',pos=111, ax_input=IR_ax, plot_pb=True)
nonIR_ax, im = make_density_map(npix,nth_neighbor,savepath=None,sources='nonIR',pos=122, ax_input=nonIR_ax, plot_pb=True)
all_ax, im = make_density_map(npix,nth_neighbor,savepath=None,sources='all_contour', pos=133, ax_input=all_ax)


IR_ax.tick_params(axis='both', direction='in', color='white')
nonIR_ax.tick_params(axis='both', direction='in', color='white')
all_ax.tick_params(axis='both', direction='in', color='white')

nonIR_ax.tick_params(axis='y', labelleft=False)
all_ax.tick_params(axis='y', labelleft=False)

IR_ax.set_ylabel('Dec.')
IR_ax.set_xlabel('RA')
nonIR_ax.set_xlabel('RA')
all_ax.set_xlabel('RA')

plt.subplots_adjust(wspace=0.07, left=0.1, right=0.9)
cbar_ax = fig.add_axes([0.91, 0.1, 0.018, 0.78]) #add_subplot(gs[3])
fig.colorbar(im, cax=cbar_ax)
cbar_ax.set_xlabel('arcsec')
#plt.tight_layout()
plt.savefig(savepath, bbox_inches='tight')



