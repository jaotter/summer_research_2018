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




def make_density_map(npix, n_neighbor, savepath):

    tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_jun20.fits')
    b3_fl = fits.open('/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
    b3_wcs = WCS(b3_fl[0].header).celestial
    b3_data = b3_fl[0].data.squeeze()
    
    b3_xmax = len(b3_data)
    b3_ymax = len(b3_data[0])

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
    nth_dists = (nth_dists_pix * pixel_scale).value
    
    density_map = np.reshape(nth_dists, (npix,npix))
    
    #make plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=new_wcs)
    im = ax.imshow(density_map, origin='lower', cmap=plt.cm.inferno_r)
    plt.colorbar(im)

    xlim = ax.get_xlim()
    ax.scatter(src_map_coords[1], src_map_coords[0], marker='*')
    ax.set_xlim(xlim)
    
    plt.savefig(savepath, bbox_inches='tight')


npix=512
nth_neighbor=3
savepath = f'/home/jotter/nrao/plots/density_map_npix{npix}.png'
make_density_map(npix,nth_neighbor,savepath)
