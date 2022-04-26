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
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy.stats

from catalogs import coupviz, forbviz, eisner2016_mm, eisner2016_gems, MLLA, meingastsci
from catalogs import forbcoords, coupcoords, mm_coords, MLLAcoords, full_mm_sourcelist, full_irxradmm_sourcelist, full_xradmm_sourcelist, IR_coords

distance = 400*u.pc

figsize=(9,10)

basepath = '/Users/adam/work/students/JustinOtter/summer_research_2018/'

npix = 400

BL_frame = SkyCoord(ra='5:35:18', dec='-5:23:30', unit=(u.hourangle,u.degree))
TR_frame = SkyCoord(ra='5:35:12', dec='-5:21:50', unit=(u.hourangle,u.degree))

for zoom, cdelt, bwfactor in (("", 0.6, 1), ("_zoom", 0.06, 0.5)):


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

    mywcs = new_wcs = WCS(new_header).celestial

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

    def rescale_densitymap(density_map, cat, surf=False):
        peak = np.unravel_index(np.argmax(density_map), xx.shape)
        maxpos = new_wcs.pixel_to_world(xx[peak], yy[peak])
        idx, sep, _ = maxpos.match_to_catalog_sky(cat, nthneighbor=nthneighbor)
        if surf:
            density_map = (nthneighbor / (np.pi * (sep*distance).to(u.pc, u.dimensionless_angles())**2)) * 0.5 * u.M_sun * density_map/density_map.max()
        else:
            density_map = (nthneighbor / (4/3*np.pi * (sep*distance).to(u.pc, u.dimensionless_angles())**3)).value * density_map / density_map.max()
        return density_map

    nthneighbor = 5

    kde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(full_mm_sourcelist), bw_method=0.2*bwfactor)
    density_map = np.reshape(kde(positions), xx.shape)
    density_map = rescale_densitymap(density_map, full_mm_sourcelist)

    coupkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(coupcoords), bw_method=0.05*bwfactor)
    coup_density_map = np.reshape(coupkde(positions), xx.shape)
    coup_density_map = rescale_densitymap(coup_density_map, coupcoords)

    forbkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(forbcoords), bw_method=0.1*bwfactor)
    forb_density_map = np.reshape(forbkde(positions), xx.shape)
    forb_density_map = rescale_densitymap(forb_density_map, forbcoords)

    jointkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(full_xradmm_sourcelist),
                                        bw_method=0.05*bwfactor)
    jointdensity_map = np.reshape(jointkde(positions), xx.shape)
    jointdensity_map = rescale_densitymap(jointdensity_map, full_xradmm_sourcelist)

    IRkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(IR_coords), bw_method=0.05*bwfactor)
    IR_density_map = np.reshape(IRkde(positions), xx.shape)
    IR_density_map = rescale_densitymap(IR_density_map, IR_coords)


    jointIRkde = scipy.stats.gaussian_kde(new_wcs.world_to_pixel(full_irxradmm_sourcelist),
                                        bw_method=0.05*bwfactor)
    jointIRdensity_map = np.reshape(jointIRkde(positions), xx.shape)
    jointIRsurfdensity_map = rescale_densitymap(jointIRdensity_map, full_irxradmm_sourcelist, surf=True)
    jointIRdensity_map = rescale_densitymap(jointIRdensity_map, full_irxradmm_sourcelist)


    fig = pl.figure(num=0, figsize=figsize)
    fig.clf()
    ax = fig.add_subplot(projection=new_wcs)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm:ss.s')
    ra.ticklabels.set_fontsize(18)
    dec.ticklabels.set_fontsize(18)

    im = ax.imshow(jointIRsurfdensity_map,
                   norm=visualization.simple_norm(jointIRsurfdensity_map, stretch='linear'),
                   cmap='gray')

    if zoom:
        B7_coord = SkyCoord(83.8104626, -5.37515542, unit=u.degree)
        Mb7 = regions.CircleSkyRegion(B7_coord, radius=10*u.arcsec)
        Mb7p = Mb7.to_pixel(ax.wcs)
        ax.axis([Mb7p.center.x - Mb7p.radius*1.1, Mb7p.center.x + Mb7p.radius*1.1,
                 Mb7p.center.y - Mb7p.radius*1.1, Mb7p.center.y + Mb7p.radius*1.1,])
    else:
        BL_pix = mywcs.all_world2pix(BL_frame.ra.degree, BL_frame.dec.degree, 0)
        TR_pix = mywcs.all_world2pix(TR_frame.ra.degree, TR_frame.dec.degree, 0)

        #ax.set_ylim(0,1250)
        #ax.set_xlim(100,1556)
        ax.set_ylim(BL_pix[1], TR_pix[1])
        ax.set_xlim(BL_pix[0], TR_pix[0])

        ax.set_xlabel('Right Ascension', fontsize=18)
        ax.set_ylabel('Declination', fontsize=18)



    divider = make_axes_locatable(ax)
    cax1 = fig.add_axes([ax.get_position().x1-0.05,
                         ax.get_position().y0,
                         0.03,
                         ax.get_position().height])
    cb = pl.colorbar(mappable=im, cax=cax1)
    cb.set_label("M$_\odot$ / pc$^{2}$", fontsize=18)
    cb.ax.tick_params(labelsize=18)
    ax.set_xlabel('Right Ascension', fontsize=18)
    ax.set_ylabel('Declination', fontsize=18)
    fig.savefig(f'{basepath}/figures/surface_density_map_joint_0.5msun.png')



    fig = pl.figure(num=1, figsize=figsize)
    fig.clf()
    ax = fig.add_subplot(projection=new_wcs)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm:ss.s')
    ra.ticklabels.set_fontsize(18)
    dec.ticklabels.set_fontsize(18)

    im = ax.imshow(density_map,
                   norm=visualization.simple_norm(density_map, stretch='linear'),
                   cmap='gray')

    if zoom:
        B7_coord = SkyCoord(83.8104626, -5.37515542, unit=u.degree)
        Mb7 = regions.CircleSkyRegion(B7_coord, radius=10*u.arcsec)
        Mb7p = Mb7.to_pixel(ax.wcs)
        ax.axis([Mb7p.center.x - Mb7p.radius*1.1, Mb7p.center.x + Mb7p.radius*1.1,
                 Mb7p.center.y - Mb7p.radius*1.1, Mb7p.center.y + Mb7p.radius*1.1,])
    else:
        BL_pix = mywcs.all_world2pix(BL_frame.ra.degree, BL_frame.dec.degree, 0)
        TR_pix = mywcs.all_world2pix(TR_frame.ra.degree, TR_frame.dec.degree, 0)

        #ax.set_ylim(0,1250)
        #ax.set_xlim(100,1556)
        ax.set_ylim(BL_pix[1], TR_pix[1])
        ax.set_xlim(BL_pix[0], TR_pix[0])

        ax.set_xlabel('Right Ascension', fontsize=18)
        ax.set_ylabel('Declination', fontsize=18)



    divider = make_axes_locatable(ax)
    cax1 = fig.add_axes([ax.get_position().x1-0.05,
                         ax.get_position().y0,
                         0.03,
                         ax.get_position().height])
    cb = pl.colorbar(mappable=im, cax=cax1)
    cb.set_label("Stars / pc$^{3}$", fontsize=18)
    cb.ax.tick_params(labelsize=18)
    ax.set_xlabel('Right Ascension', fontsize=18)
    ax.set_ylabel('Declination', fontsize=18)
    pl.savefig(f'{basepath}/figures/density_map_mm{zoom}.png', )

    #ax.contour(density_map, cmap='gray')
    cn = ax.contour(forb_density_map, cmap='inferno')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))

    pl.savefig(f'{basepath}/figures/density_map_mm_radio{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)

    cn = ax.contour(coup_density_map, cmap='viridis')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_mm_xray{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)



    cn = ax.contour(IR_density_map, cmap='magma')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_mm_IR{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)





    #ax.contour(jointdensity_map, cmap='jet')
    fig = pl.figure(num=2, figsize=figsize)
    fig.clf()
    ax = fig.add_subplot(projection=new_wcs)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm:ss.s')
    ra.ticklabels.set_fontsize(18)
    dec.ticklabels.set_fontsize(18)


    im = ax.imshow(jointdensity_map,
                   norm=visualization.simple_norm(jointdensity_map, stretch='linear'),
                   cmap='gray')

    if zoom:
        B7_coord = SkyCoord(83.8104626, -5.37515542, unit=u.degree)
        Mb7 = regions.CircleSkyRegion(B7_coord, radius=10*u.arcsec)
        Mb7p = Mb7.to_pixel(ax.wcs)
        ax.axis([Mb7p.center.x - Mb7p.radius*1.1, Mb7p.center.x + Mb7p.radius*1.1,
                 Mb7p.center.y - Mb7p.radius*1.1, Mb7p.center.y + Mb7p.radius*1.1,])
    else:
        BL_pix = mywcs.all_world2pix(BL_frame.ra.degree, BL_frame.dec.degree, 0)
        TR_pix = mywcs.all_world2pix(TR_frame.ra.degree, TR_frame.dec.degree, 0)

        #ax.set_ylim(0,1250)
        #ax.set_xlim(100,1556)
        ax.set_ylim(BL_pix[1], TR_pix[1])
        ax.set_xlim(BL_pix[0], TR_pix[0])

        ax.set_xlabel('Right Ascension', fontsize=18)
        ax.set_ylabel('Declination', fontsize=18)

    divider = make_axes_locatable(ax)
    cax1 = fig.add_axes([ax.get_position().x1-0.05,
                         ax.get_position().y0,
                         0.03,
                         ax.get_position().height])
    cb = pl.colorbar(mappable=im, cax=cax1)
    cb.set_label("Stars / pc$^{3}$", fontsize=18)
    cb.ax.tick_params(labelsize=18)
    ax.set_xlabel('Right Ascension', fontsize=18)
    ax.set_ylabel('Declination', fontsize=18)
    pl.savefig(f'{basepath}/figures/density_map_merged{zoom}.png', )

    cn = ax.contour(forb_density_map, cmap='inferno')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))

    pl.savefig(f'{basepath}/figures/density_map_merged_radio{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)

    cn = ax.contour(coup_density_map, cmap='viridis')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_merged_xray{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)

    cn = ax.contour(density_map, cmap='jet')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_merged_mm{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)


    cn = ax.contour(IR_density_map, cmap='magma')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_merged_IR{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)





    fig = pl.figure(num=3, figsize=figsize)
    fig.clf()
    ax = fig.add_subplot(projection=new_wcs)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm:ss.s')
    ra.ticklabels.set_fontsize(18)
    dec.ticklabels.set_fontsize(18)


    im = ax.imshow(jointIRdensity_map,
                   norm=visualization.simple_norm(jointIRdensity_map, stretch='linear'),
                   cmap='gray')

    if zoom:
        B7_coord = SkyCoord(83.8104626, -5.37515542, unit=u.degree)
        Mb7 = regions.CircleSkyRegion(B7_coord, radius=10*u.arcsec)
        Mb7p = Mb7.to_pixel(ax.wcs)
        ax.axis([Mb7p.center.x - Mb7p.radius*1.1, Mb7p.center.x + Mb7p.radius*1.1,
                 Mb7p.center.y - Mb7p.radius*1.1, Mb7p.center.y + Mb7p.radius*1.1,])
    else:
        BL_pix = mywcs.all_world2pix(BL_frame.ra.degree, BL_frame.dec.degree, 0)
        TR_pix = mywcs.all_world2pix(TR_frame.ra.degree, TR_frame.dec.degree, 0)

        #ax.set_ylim(0,1250)
        #ax.set_xlim(100,1556)
        ax.set_ylim(BL_pix[1], TR_pix[1])
        ax.set_xlim(BL_pix[0], TR_pix[0])

        ax.set_xlabel('Right Ascension', fontsize=18)
        ax.set_ylabel('Declination', fontsize=18)

    divider = make_axes_locatable(ax)
    cax1 = fig.add_axes([ax.get_position().x1-0.05,
                         ax.get_position().y0,
                         0.03,
                         ax.get_position().height])
    cb = pl.colorbar(mappable=im, cax=cax1)
    cb.set_label("Stars / pc$^{3}$", fontsize=18)
    cb.ax.tick_params(labelsize=18)
    ax.set_xlabel('Right Ascension', fontsize=18)
    ax.set_ylabel('Declination', fontsize=18)
    pl.savefig(f'{basepath}/figures/density_map_mergedIR{zoom}.png', )

    cn = ax.contour(forb_density_map, cmap='inferno')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))

    pl.savefig(f'{basepath}/figures/density_map_mergedIR_radio{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)

    cn = ax.contour(coup_density_map, cmap='viridis')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_mergedIR_xray{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)

    cn = ax.contour(density_map, cmap='jet')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_mergedIR_mm{zoom}.png', )
    for cll in cn.collections + axhls:
        cll.set_visible(False)


    cn = ax.contour(IR_density_map, cmap='magma')
    cax = cb.ax
    axhls = []
    for level, coll in zip(cn.levels, cn.collections):
        axhls.append(cax.axhline(level, color=coll.get_color()[0]))
    pl.savefig(f'{basepath}/figures/density_map_mergedIR_IR{zoom}.png', )
    #for cll in cn.collections + axhls:
    #    cll.set_visible(False)


    pl.savefig(f'{basepath}/figures/density_map_mergedIR_IR_zoom{zoom}.png', )
