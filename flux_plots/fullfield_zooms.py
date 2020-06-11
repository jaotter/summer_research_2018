import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
import astropy.visualization
from astropy.convolution import convolve, Gaussian2DKernel
from mpl_plot_templates import asinh_norm
import matplotlib
from collections import defaultdict
import warnings
import astropy.units as u
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


def calc_bbox(center, sidelength):
    BL_ra = center.ra.degree + sidelength/2
    BL_dec = center.dec.degree - sidelength/2
    TR_ra = center.ra.degree - sidelength/2
    TR_dec = center.dec.degree + sidelength/2
    BL = coordinates.SkyCoord(BL_ra, BL_dec, unit=u.degree)
    TR = coordinates.SkyCoord(TR_ra, TR_dec, unit=u.degree)
    return BL, TR

def create_zoomregion(srcid, table, bbox, locs=[1,3], name=None, vmin=-0.0001):
    srcind = np.where(table['D_ID'] == srcid)[0]
    srctab = table[srcind]
    src_coord = coordinates.SkyCoord(ra=srctab['RA_B3'].data, dec=srctab['DEC_B3'].data, unit=u.degree)
    src_fwhm = srctab['fwhm_maj_B3']*u.arcsecond
    sidelength = src_fwhm.to(u.degree)*4
    bottomleft, topright = calc_bbox(src_coord, sidelength.value)

    zoom = 10
    vmax = 90

    #keys in reg_dict: 'bottomleft':SkyCoord, 'topright':SkyCoord, 'bbox':[float, float] (location of inset BL),
    #l1,l2: int - corner to draw line (UR is 1, then ccw), min:float - vmin, max:float - vmax, if greater than 1 than is percent of max flux in inset, zoom:float
    
    key = f'{name+" " if name is not None else ""}({srcid})'
    reg_dict = {'bottomleft':bottomleft, 'topright':topright, 'bbox':bbox, 'min':vmin, 'max':vmax, 'zoom':zoom, 'loc':2, 'l1':locs[0], 'l2':locs[1]}

    return key, reg_dict


table = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_apr20.fits')

#INSET 1
#sources = [30, 43, 31, 23, 38, 32, 28, 45, 20, 71, 22, 36]
#bboxes = [[0.25, 0.9],[0.725,0.85],[0.25,0.5],[0.75,0.25],[0.75,0.65],[0.75,0.50],[0.4,0.27],[0.55,0.85],[0.53,0.27],[0.65,0.25],[0.25,0.3],[0.25,0.7]]
#names = ['Source I', 'BN', None,'Source N','IRC6E','IRC2C',None,None,None,None,None,None]
#locs_all = [[1,3],[2,3],[2,4],[1,3],[2,3],[1,3],[1,2],[1,4],[1,2],[1,2],[2,4],[1,4]]
#filename = 'B3_inset1.png'
#vmin = -0.0001
    
#INSET 2
sources = [34, 39, 41, 76, 79, 29, 80, 75, 35, 42, 49, 50]
bboxes = [[0.6,0.4],[0.75,0.7],[0.5,0.72],[0.5,0.4],[0.6,0.55],[0.25,0.5],[0.4,0.27],[0.3,0.3],[0.25,0.66],[0.25, 0.9],[0.45,0.9],[0.65,0.9]]
names = np.repeat(None, len(sources))
locs_all = [[1,3],[2,3],[2,3],[1,3],[2,3],[1,4],[1,2],[2,4],[1,3],[3,4],[2,4],[2,3]]
filename = 'B3_inset2.png'
vmin = -0.00001
    
zoomregions_auto = {}

for i in range(len(sources)):
    key, reg_dict = create_zoomregion(sources[i], table, bboxes[i], locs=locs_all[i], name=names[i], vmin=vmin)
    zoomregions_auto[key] = reg_dict

srcI_ind = np.where(table['D_ID'] == 30)[0]
srcI_coord = coordinates.SkyCoord(table['RA_B3'][srcI_ind], table['DEC_B3'][srcI_ind], unit=u.degree)
BL, TR = calc_bbox(srcI_coord, sidelength=(30*u.arcsec).to(u.degree).value)

#zoomregions_auto['(72,37)'] = {'bottomleft': coordinates.SkyCoord(ra=["5:35:14.435"],
#                                                   dec=["-5:22:28.55"],
#                                                   unit=(u.h, u.deg),
#                                                   frame='icrs'),
#                'topright': coordinates.SkyCoord(ra=["5:35:14.403"],
#                                                 dec=["-5:22:28.28"],
#                                                 unit=(u.h, u.deg),
#                                                 frame='icrs'),
#                'bbox':[0.4,0.85],'loc': 2,'l1':3,'l2':1,'min': -0.0001,'max': 0.0005,'zoom': 10}


zoomregions = zoomregions_auto

#psf_center = coordinates.SkyCoord('5:35:14.511 -5:22:30.561',
#                                  unit=(u.h, u.deg),
#                                  frame='icrs')


def inset_overlays(fn, zoomregions, fignum=1,
                   psffn=None,
                   vmin=-0.001, vmax=0.01,
                   directory = '/lustre/aoc/students/jotter/directory/',
                   bottomleft=coordinates.SkyCoord('5:35:15.236', '-5:22:41.85', unit=(u.h, u.deg), frame='icrs'),
                   topright=coordinates.SkyCoord('5:35:13.586', '-5:22:17.12', unit=(u.h, u.deg), frame='icrs'),
                   tick_fontsize=pl.rcParams['axes.labelsize']):
                  
    fn = directory+fn
    hdu = fits.open(fn)[0]

    mywcs = wcs.WCS(hdu.header).celestial

    figure = pl.figure(fignum)
    figure.clf()
    ax = figure.add_axes([0.15, 0.1, 0.8, 0.8], projection=mywcs)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    ra.set_axislabel("RA (ICRS)", fontsize=pl.rcParams['axes.labelsize'])
    dec.set_axislabel("Dec (ICRS)", fontsize=pl.rcParams['axes.labelsize'], minpad=0.0)
    ra.ticklabels.set_fontsize(tick_fontsize)
    ra.set_ticks(exclude_overlapping=True)
    dec.ticklabels.set_fontsize(tick_fontsize)
    dec.set_ticks(exclude_overlapping=True)


    im = ax.imshow(hdu.data.squeeze()*1e3,
                   transform=ax.get_transform(mywcs),
                   vmin=vmin*1e3, vmax=vmax*1e3, cmap=pl.cm.gray_r,
                   interpolation='nearest',
                   origin='lower', norm=asinh_norm.AsinhNorm())
    
    (x1,y1),(x2,y2) = (mywcs.wcs_world2pix([[bottomleft.ra.deg[0],
                                             bottomleft.dec.deg[0]]],0)[0],
                       mywcs.wcs_world2pix([[topright.ra.deg[0],
                                             topright.dec.deg[0]]],0)[0]
                      )

    # we'll want this later
    #make_scalebar(ax, scalebarpos,
    #              length=(0.5*u.pc / distance).to(u.arcsec,
    #                                              u.dimensionless_angles()),
    #              color='k',
    #              label='0.5 pc',
    #              text_offset=1.0*u.arcsec,
    #             )


    ax.set_aspect(1)
    ax.axis([x1,x2,y1,y2])


    for zoomregion in zoomregions:

        ZR = zoomregions[zoomregion]

        parent_ax = zoomregions[ZR['inset_axes']]['axins'] if 'inset_axes' in ZR else ax

        bl, tr = ZR['bottomleft'],ZR['topright'],
        
        (zx1,zy1),(zx2,zy2) = (mywcs.wcs_world2pix([[bl.ra.deg[0],
                                                     bl.dec.deg[0]]],0)[0],
                               mywcs.wcs_world2pix([[tr.ra.deg[0],
                                                     tr.dec.deg[0]]],0)[0]
                              )
        print(zoomregion,zx1,zy1,zx2,zy2)

        inset_data = hdu.data.squeeze()[int(zy1):int(zy2), int(zx1):int(zx2)]
        #inset_data = hdu.data.squeeze()
        inset_wcs = mywcs.celestial[int(zy1):int(zy2), int(zx1):int(zx2)]
        #inset_wcs = mywcs

        axins = zoomed_inset_axes(parent_ax, zoom=ZR['zoom'], loc=ZR['loc'],
                                  bbox_to_anchor=ZR['bbox'],
                                  bbox_transform=figure.transFigure,
                                  axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                  axes_kwargs=dict(wcs=inset_wcs))
        ZR['axins'] = axins

        vmax = ZR['max']
        if vmax > 1:
            vmax = (vmax/100)*np.nanmax(inset_data)
            
        imz = axins.imshow(inset_data,
                           #transform=parent_ax.get_transform(inset_wcs),
                           extent=[int(zx1), int(zx2), int(zy1), int(zy2)],
                           vmin=ZR['min'], vmax=vmax, cmap=pl.cm.gray_r,
                           interpolation='nearest',
                           origin='lower', norm=asinh_norm.AsinhNorm())


        ax.axis([x1,x2,y1,y2])
        #axins.axis([zx1,zx2,zy1,zy2])
        #print(axins.axis())

        axins.set_xticklabels([])
        axins.set_yticklabels([])

        #parent_ax.text(zx1, zy1, zoomregion)
        axins.text(int(zx1)+5, int(zy1)+5, zoomregion, fontsize=6)
        lon = axins.coords['ra']
        lat = axins.coords['dec']
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)

        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(parent_axes=parent_ax, inset_axes=axins,
                   loc1=ZR['l1'], loc2=ZR['l2'], fc="none", ec="0.5",
                   lw=0.5)


        figure.canvas.draw()
        assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

    cax = figure.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                           ax.bbox._bbox.y1-ax.bbox._bbox.y0])
    cb = figure.colorbar(mappable=im, cax=cax)
    #print("1. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.set_label("$S_{1 mm}$ [mJy beam$^{-1}$]")
    #print("2. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.formatter.format = "%3.1f"
    #print("3. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.set_ticks(cb.formatter.locs)
    #print("4. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.set_ticklabels(["{0:3.1f}".format(float(x)) for x in cb.formatter.locs])
    #print("5. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.ax.set_yticklabels(["{0:3.1f}".format(float(x.get_text())) for x in cb.ax.get_yticklabels()])
    #print("6. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    


    if psffn is not None:
        psf = fits.open(psffn)
        psfwcs = wcs.WCS(psf[0].header)
        cx,cy = psfwcs.celestial.wcs_world2pix(psf_center.ra.deg, psf_center.dec.deg, 0)
        cx = int(cx)
        cy = int(cy)
        zy1 = cy-50
        zy2 = cy+50
        zx1 = cx-50
        zx2 = cx+50

        inset_wcs = psfwcs.celestial[zy1:zy2, zx1:zx2]
        inset_data = psf[0].data[cy-50:cy+50, cx-50:cx+50]

        axins = zoomed_inset_axes(parent_ax, zoom=10, loc=2,
                                  bbox_to_anchor=(0.05,0.25),
                                  bbox_transform=figure.transFigure,
                                  axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                  axes_kwargs=dict(wcs=inset_wcs),
                                 )
        imz = axins.imshow(inset_data,
                           extent=[int(zx1), int(zx2), int(zy1), int(zy2)],
                           vmin=0, vmax=1, cmap=pl.cm.gray_r,
                           interpolation='nearest',
                           origin='lower', norm=asinh_norm.AsinhNorm())
        axins.contour(np.linspace(zx1, zx2, inset_data.shape[1]),
                      np.linspace(zy1, zy2, inset_data.shape[0]),
                      inset_data,
                      levels=[0.05, 0.1, 0.2, 0.3],
                      linewidths=[0.3]*10,
                      alpha=0.75,
                      #colors=['r']*10,
                     )
        axins.set_xticks([])
        axins.set_yticks([])
        axins.set_xticklabels([])
        axins.set_yticklabels([])
        lon = axins.coords['ra']
        lat = axins.coords['dec']
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)


    return figure


if __name__ == "__main__":
    fn = 'Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'
    directory='/home/jotter/nrao/images/'
    figure = inset_overlays(fn, zoomregions=zoomregions, directory=directory,
                                vmin=-0.0005, vmax=0.003, bottomleft=BL, topright=TR)
    figure.savefig('/home/jotter/nrao/plots/inset_plots/'+filename, bbox_inches='tight', dpi=300)
