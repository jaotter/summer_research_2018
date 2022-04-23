# new version, following Forbrich...
import os
from astroquery.vizier import Vizier
import regions
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import coordinates
from astropy import visualization
from astropy.wcs import WCS
from scipy.spatial import KDTree
from matplotlib.patches import Circle
import astropy.wcs.utils as wcs_utils
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import scipy.stats

distance = 400*u.pc

basepath = '/Users/adam/work/students/JustinOtter/summer_research_2018/'

def getcat(catname):
    filename = os.path.join(basepath, 'tables', catname.replace("/","_")+".ecsv")
    if os.path.exists(filename):
        return Table.read(filename)
    else:
        result = Vizier(row_limit=100000).get_catalogs(catname)[0]
        result.write(filename)
        return result

#coupviz = Vizier(row_limit=100000).get_catalogs('J/ApJS/160/319/coup')[0]
#forbviz = Vizier(row_limit=100000).get_catalogs('J/ApJ/822/93')[0]
#eisner2016_gems = Vizier(row_limit=100000).get_catalogs('J/ApJ/826/16/table1')[0]
#eisner2016_mm = Vizier(row_limit=100000).get_catalogs('J/ApJ/826/16/table2')[0]
#MLLA = Vizier(row_limit=100000).get_catalogs('J/ApJ/573/366')[0]
coupviz = getcat('J/ApJS/160/319/coup')
forbviz = getcat('J/ApJ/822/93')
eisner2016_gems = getcat('J/ApJ/826/16/table1')
eisner2016_mm = getcat('J/ApJ/826/16/table2')
MLLA = getcat('J/ApJ/573/366')

forbcoords = SkyCoord([x for x in forbviz['RAJ2000']],
                      [x for x in forbviz['DEJ2000']], frame='fk5', unit=(u.h,u.deg))
forbregs = regions.Regions([regions.PointSkyRegion(coord, meta={'text': id, 'shape': 'square', 'color': 'orange'})
                            for coord,id in zip(forbcoords, forbviz['Seq'])])
forbregs.write(f'{basepath}/tables/Forbrich2016.reg', overwrite=True)

coupcoords = SkyCoord([x for x in coupviz['RAJ2000']],
                      [x for x in coupviz['DEJ2000']], frame='fk5', unit=(u.deg, u.deg))

coupregs = regions.Regions([regions.PointSkyRegion(coord, meta={'text': id, 'shape': 'diamond', 'color': 'cyan'})
                            for coord,id in zip(coupcoords, coupviz['COUP'])])
coupregs.write(f'{basepath}/tables/COUP.reg', overwrite=True)


eisner2018_mm = Table.read(f'{basepath}/tables/eisner_tbl.txt', format='ascii')

crd_eis2016mm = SkyCoord(eisner2016_mm['RAJ2000'], eisner2016_mm['DEJ2000'], frame='fk5', unit=(u.h, u.deg))
crd_eis2016gems = SkyCoord(eisner2016_gems['RAJ2000'], eisner2016_gems['DEJ2000'], frame='fk5', unit=(u.h, u.deg))
crd_eis2018mm = SkyCoord(eisner2018_mm['alpha'], eisner2018_mm['delta'], frame='fk5', unit=(u.h, u.deg))

tab = B3_table = Table.read(f'{basepath}/final_tables/datafile4.txt', format='ascii.cds')
mm_coords = SkyCoord((B3_table['RAs'].quantity + (5*u.hour).to(u.s) + (35*u.min).to(u.s))*(15*u.arcsec/u.s),
                  -5*u.deg - B3_table['DEm'].quantity - B3_table['DEs'].quantity, frame='icrs', )

idx, sep, _ = crd_eis2018mm.match_to_catalog_sky(mm_coords)
new_in_e2018 = sep > 1*u.arcsec

full_mm_sourcelist = coordinates.concatenate([mm_coords.fk5, crd_eis2018mm[new_in_e2018]])
full_mm_sourcelist = SkyCoord(full_mm_sourcelist.ra, full_mm_sourcelist.dec, frame='fk5')

idx, sep, _ = crd_eis2016mm.match_to_catalog_sky(full_mm_sourcelist)
new_in_e2016 = sep > 1*u.arcsec

full_mm_sourcelist = coordinates.concatenate([full_mm_sourcelist.fk5,
                                              crd_eis2016mm[new_in_e2016].fk5])
full_mm_sourcelist = SkyCoord(full_mm_sourcelist.ra, full_mm_sourcelist.dec, frame='fk5')



idx, sep, _ = forbcoords.match_to_catalog_sky(full_mm_sourcelist)
new_in_forb = sep > 1*u.arcsec
full_radmm_sourcelist = coordinates.concatenate([full_mm_sourcelist.fk5,
                                                 forbcoords[new_in_forb].fk5])
full_radmm_sourcelist = SkyCoord(full_radmm_sourcelist.ra, full_radmm_sourcelist.dec, frame='fk5')



idx, sep, _ = coupcoords.match_to_catalog_sky(full_mm_sourcelist)
new_in_coup = sep > 1*u.arcsec
full_xradmm_sourcelist = coordinates.concatenate([full_radmm_sourcelist.fk5,
                                                  coupcoords[new_in_coup].fk5])
full_xradmm_sourcelist = SkyCoord(full_xradmm_sourcelist.ra, full_xradmm_sourcelist.dec, frame='fk5')



MLLAcoords = SkyCoord(MLLA['RAJ2000'], MLLA['RAJ2000'], frame='fk5', unit=(u.h, u.deg))

idx, sep, _ = crd_eis2016gems.match_to_catalog_sky(MLLAcoords)
new_in_e2016gems = sep > 1*u.arcsec
IR_coords = coordinates.concatenate([MLLAcoords.fk5,
                                     crd_eis2016gems[new_in_e2016gems].fk5])
IR_coords = SkyCoord(IR_coords.ra, IR_coords.dec, frame='fk5')





npix = 400
cdelt = npix / 600



#creating new wcs for density map
new_header = {}
new_header['CDELT1'] = -cdelt/3600
new_header['CDELT2'] = cdelt/3600
new_header['CRPIX1'] = npix/2
new_header['CRPIX2'] = npix/2
new_header['CRVAL1'] = 83.8104973
new_header['CRVAL2'] = -5.3751787
new_header['NAXIS1'] = npix
new_header['NAXIS2'] = npix
new_header['NAXIS'] = 2
new_header['CTYPE1'] = 'RA---TAN'
new_header['CTYPE2'] = 'DEC--TAN'

new_wcs = WCS(new_header).celestial

#density_map = np.zeros((npix,npix))
xx, yy = np.mgrid[0:npix, 0:npix] #grid with density map pixel coords
positions = np.vstack([xx.ravel(), yy.ravel()])

pix_points = np.array(list(zip(xx.ravel()+0.5,yy.ravel()+0.5)))

mm_pix_coords = np.array(new_wcs.world_to_pixel(full_mm_sourcelist))
mm_map_coords = mm_pix_coords
mm_pix_coords = np.transpose(mm_pix_coords)

"""
tree = KDTree(src_pix_coords) #KDTree has source positions in map pix coords

n_neighbor = 10
full_dists, full_inds = tree.query(pix_points, k=n_neighbor) #each entry of full_dists has nearest neighbors up to the nth nearest

nth_dists_pix = [dist[n_neighbor-1] for dist in full_dists]
"""

#convert to arcseconds
pixel_scales = wcs_utils.proj_plane_pixel_scales(new_wcs)

pixel_scale = (pixel_scales[0]*u.degree).to(u.arcsecond)
#nth_dists = (nth_dists_pix * pixel_scale).value

#density_map = np.reshape(1./nth_dists, (npix,npix))

def rescale_densitymap(density_map, cat):
    peak = np.unravel_index(np.argmax(density_map), xx.shape)
    maxpos = new_wcs.pixel_to_world(xx[peak], yy[peak])
    idx, sep, _ = maxpos.match_to_catalog_sky(cat, nthneighbor=nthneighbor)
    density_map = (nthneighbor / (4/3*np.pi * (sep*distance).to(u.pc, u.dimensionless_angles())**3)).value * density_map / density_map.max()
    return density_map

nthneighbor = 5

kde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(full_mm_sourcelist), bw_method=0.2)
density_map = np.reshape(kde(positions), xx.shape)
density_map = rescale_densitymap(density_map, full_mm_sourcelist)

coupkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(coupcoords), bw_method=0.05)
coup_density_map = np.reshape(coupkde(positions), xx.shape)
coup_density_map = rescale_densitymap(coup_density_map, coupcoords)

forbkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(forbcoords), bw_method=0.1)
forb_density_map = np.reshape(forbkde(positions), xx.shape)
forb_density_map = rescale_densitymap(forb_density_map, forbcoords)

jointkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(full_xradmm_sourcelist),
                                    bw_method=0.05)
jointdensity_map = np.reshape(jointkde(positions), xx.shape)
peak = np.unravel_index(np.argmax(jointdensity_map), xx.shape)
jointdensity_map = rescale_densitymap(jointdensity_map, full_xradmm_sourcelist)

IRkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(IR_coords), bw_method=0.05)
IR_density_map = np.reshape(IRkde(positions), xx.shape)
IR_density_map = rescale_densitymap(IR_density_map, IR_coords)




header = fits.Header.fromtextfile(f'{basepath}/trapezium_small_photoshop.wcs')
mywcs = WCS(header).celestial


fig = pl.figure(num=1, figsize=(10,10))
fig.clf()
ax = fig.add_subplot(projection=new_wcs)
im = ax.imshow(density_map,
               norm=visualization.simple_norm(density_map, stretch='linear'),
               cmap='gray')
cb = pl.colorbar(mappable=im)
cb.set_label("Stars / pc$^{3}$")
pl.savefig(f'{basepath}/figures/density_map_mm.png', bbox_inches='tight')
ax.set_xlabel('Right Ascension', fontsize=16)
ax.set_ylabel('Declination', fontsize=16)

#ax.contour(density_map, cmap='gray')
cn = ax.contour(forb_density_map, cmap='inferno')
cax = cb.ax
axhls = []
for level, coll in zip(cn.levels, cn.collections):
    axhls.append(cax.axhline(level, color=coll.get_color()[0]))

pl.savefig(f'{basepath}/figures/density_map_mm_radio.png', bbox_inches='tight')
for cll in cn.collections + axhls:
    cll.set_visible(False)

cn = ax.contour(coup_density_map, cmap='viridis')
cax = cb.ax
axhls = []
for level, coll in zip(cn.levels, cn.collections):
    axhls.append(cax.axhline(level, color=coll.get_color()[0]))
pl.savefig(f'{basepath}/figures/density_map_mm_xray.png', bbox_inches='tight')
for cll in cn.collections + axhls:
    cll.set_visible(False)

#ax.contour(jointdensity_map, cmap='jet')
fig = pl.figure(num=2, figsize=(10,10))
fig.clf()
ax = fig.add_subplot(projection=new_wcs)
im = ax.imshow(jointdensity_map,
               norm=visualization.simple_norm(jointdensity_map, stretch='linear'),
               cmap='gray')
cb = pl.colorbar(mappable=im)
cb.set_label("Stars / pc$^{3}$")
ax.set_xlabel('Right Ascension', fontsize=16)
ax.set_ylabel('Declination', fontsize=16)
pl.savefig(f'{basepath}/figures/density_map_merged.png', bbox_inches='tight')

cn = ax.contour(forb_density_map, cmap='inferno')
cax = cb.ax
axhls = []
for level, coll in zip(cn.levels, cn.collections):
    axhls.append(cax.axhline(level, color=coll.get_color()[0]))

pl.savefig(f'{basepath}/figures/density_map_merged_radio.png', bbox_inches='tight')
for cll in cn.collections + axhls:
    cll.set_visible(False)

cn = ax.contour(coup_density_map, cmap='viridis')
cax = cb.ax
axhls = []
for level, coll in zip(cn.levels, cn.collections):
    axhls.append(cax.axhline(level, color=coll.get_color()[0]))
pl.savefig(f'{basepath}/figures/density_map_merged_xray.png', bbox_inches='tight')
for cll in cn.collections + axhls:
    cll.set_visible(False)

cn = ax.contour(density_map, cmap='jet')
cax = cb.ax
axhls = []
for level, coll in zip(cn.levels, cn.collections):
    axhls.append(cax.axhline(level, color=coll.get_color()[0]))
pl.savefig(f'{basepath}/figures/density_map_merged_mm.png', bbox_inches='tight')
for cll in cn.collections + axhls:
    cll.set_visible(False)

