import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
import astropy.visualization
from astropy.convolution import convolve, Gaussian2DKernel
#from mpl_plot_templates import asinh_norm
import matplotlib
from collections import defaultdict
import warnings
import astropy.units as u
from matplotlib.patches import RegularPolygon
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
    srcind = np.where(table['ID'] == srcid)[0]
    srctab = table[srcind]
    src_coord = coordinates.SkyCoord(ra=srctab['RA_B3'].data, dec=srctab['DEC_B3'].data, unit=u.degree)
    src_fwhm = srctab['fwhm_maj_B3']*u.arcsecond
    sidelength = src_fwhm.to(u.degree)*4
    bottomleft, topright = calc_bbox(src_coord, sidelength.value)
    zoom = 10
    vmax = 90

    #keys in reg_dict: 'bottomleft':SkyCoord, 'topright':SkyCoord, 'bbox':[float, float] (location of inset BL),
    #l1,l2: int - corner to draw line (UR is 1, then ccw), min:float - vmin, max:float - vmax, if greater than 1 than is percent of max flux in inset, zoom:float

    key = f'{name if name != "     " else ""}{" " if name != "     " else ""}({srcid})'
    reg_dict = {'bottomleft':bottomleft, 'topright':topright, 'bbox':bbox, 'min':vmin, 'max':vmax, 'zoom':zoom, 'loc':2, 'l1':locs[0], 'l2':locs[1]}

    return key, reg_dict



psf_center = coordinates.SkyCoord('5:35:14.511 -5:22:30.561',
                                  unit=(u.h, u.deg),
                                  frame='icrs')


def inset_overlays(fn, zoomregions, fignum=1,
                   savefn='',
                   psffn=None,
                   vmin=-0.001, vmax=0.01,
                   directory = '.',
                   bottomleft=coordinates.SkyCoord('5:35:15.236', '-5:22:41.85', unit=(u.h, u.deg), frame='icrs'),
                   topright=coordinates.SkyCoord('5:35:13.586', '-5:22:17.12', unit=(u.h, u.deg), frame='icrs'),
                   tick_fontsize=pl.rcParams['axes.labelsize'],
                   make_overview=False,
                   colorbar=True,
                  ):

    fn = f'{directory}/images/{fn}'
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
                   origin='lower')#, norm=asinh_norm.AsinhNorm())

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

    tab = Table.read(f'{directory}/tables/r0.5_catalog_bgfit_may21_ulim.fits')
    srcI_ind = np.where(tab['D_ID'] == 30)[0]
    BN_ind = np.where(tab['D_ID'] == 43)[0]
    srcI_coord = mywcs.wcs_world2pix([[tab['RA_B3'][srcI_ind][0], tab['DEC_B3'][srcI_ind][0]]],0)
    BN_coord = mywcs.wcs_world2pix([[tab['RA_B3'][BN_ind][0], tab['DEC_B3'][BN_ind][0]]],0)

    srcI_reg = RegularPolygon((srcI_coord[0][0], srcI_coord[0][1]), 3, color='tab:cyan', fill=False, radius=75)
    BN_reg = RegularPolygon((BN_coord[0][0], BN_coord[0][1]), 4, color='tab:orange', fill=False, radius=75)

    ax.add_patch(srcI_reg)
    ax.add_patch(BN_reg)

    if make_overview:
        lon = ax.coords[0]
        lat = ax.coords[1]
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)
        figure.savefig(f'Orion_Overview_{savefn}.png', bbox_inches='tight', pad_inches=0, dpi=300)
        sc = ax.scatter(tab['RA_B3'], tab['DEC_B3'], marker='o', facecolors='none', edgecolors='r', transform=ax.get_transform('world'))
        figure.savefig(f'Orion_Overview_withsources_{savefn}.png', bbox_inches='tight', pad_inches=0, dpi=300)
        sc.set_visible(False)


    else:
        for zoomregion in zoomregions:

            ZR = zoomregions[zoomregion]

            parent_ax = zoomregions[ZR['inset_axes']]['axins'] if 'inset_axes' in ZR else ax

            bl, tr = ZR['bottomleft'],ZR['topright'],
            (zx1,zy1),(zx2,zy2) = (mywcs.wcs_world2pix([[bl.ra.deg[0],
                                                         bl.dec.deg[0]]],0)[0],
                                   mywcs.wcs_world2pix([[tr.ra.deg[0],
                                                         tr.dec.deg[0]]],0)[0])

            pix_coords = [zx1,zy1,zx2,zy2]
            ymax = len(hdu.data.squeeze())
            xmax = len(hdu.data.squeeze()[0])
            if zx1 < 0:
                zx1 = 0
            if zy1 < 0:
                zy1 = 0
            if zy2 > ymax:
                zy2 = ymax
            if zx2 > xmax:
                zx2 = xmax

            if zy2 < 0 or zx2 < 0:
                print("FAILED on ",zoomregion,zx1,zy1,zx2,zy2)
                continue

            print(zoomregion,zx1,zy1,zx2,zy2)

            inset_data = hdu.data.squeeze()[int(zy1):int(zy2), int(zx1):int(zx2)]
            #inset_data = hdu.data.squeeze()
            inset_wcs = mywcs.celestial[int(zy1):int(zy2), int(zx1):int(zx2)]

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
                               origin='lower')#, norm=asinh_norm.AsinhNorm())


            ax.axis([x1,x2,y1,y2])
            #axins.axis([zx1,zx2,zy1,zy2])
            #print(axins.axis())

            axins.set_xticklabels([])
            axins.set_yticklabels([])

            #parent_ax.text(zx1, zy1, zoomregion)
            if 'None ' in zoomregion:
                zoomregion = zoomregion.split('None ')[-1]
            axins.text(int(zx1)+5, int(zy1)+5, zoomregion, fontsize=5)
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

    if colorbar:
        cax = figure.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                               ax.bbox._bbox.y1-ax.bbox._bbox.y0])
        cb = figure.colorbar(mappable=im, cax=cax)
        #print("1. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
        cb.set_label("log $S_{3 mm}$ [mJy beam$^{-1}$]")
        #print("2. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
        #cb.formatter.format = "%3.1f"
        #print("3. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
        #cb.set_ticks(cb.formatter.locs)
        #print("4. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
        #cb.set_ticklabels(["{0:3.1f}".format(float(x)) for x in cb.formatter.locs])
        #print("5. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
        #cb.ax.set_yticklabels(["{0:3.1f}".format(float(x.get_text())) for x in cb.ax.get_yticklabels()])
        #print("6. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    else:
        lon = ax.coords[0]
        lat = ax.coords[1]
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)
        figure.savefig(f'Orion_Overview_insets_{savefn}.png', bbox_inches='tight', pad_inches=0, dpi=300)



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


    figure.savefig(savefn+".pdf", bbox_inches='tight')

    return figure


if __name__ == "__main__":

    import socket
    if 'ufhpc' in socket.gethostname():
        root = '/orange/adamginsburg/orion/jaotter_git/'
    else:
        root='/home/jotter/nrao/summer_research_2018/'

    table = Table.read(f'{root}/tables/r0.5_catalog_bgfit_may21_ulim.fits')

    all_sources = (
        [12, 69, 15, 14, 10, 9, 16, 23, 70, 29, 68, 65, 25, 31, 38, 34, 35, 39, 40, 45, 44, 21, 28, 37, 32, 20, 26, 27, 67, 30, 66, 24],
        [ 3,  4,  6,  7,  8, 11, 13, 19, 22, 36, 43, 46, 47, 48, 49, 50, 51, 53, 56, 58, 71],
        [10, 8, 13, 84, 54, 55, 53, 46],
        [5, 6, 9, 15, 40, 51, 21, 24, 48, 52, 63, 56, 61, 58, 7, 4, 2, 3, 19, 47, 57, 60, 62, 59, 67, 66, 64, 70],
        [30, 43, 31, 23, 38, 32, 28, 45, 20, 71, 22, 36],
        [34, 39, 41, 76, 79, 29, 80, 75, 35, 42, 49, 50],
        [14, 81, 17, 16, 12, 11, 18, 27, 83, 33],
    )
    all_bboxes = (
              [[0.34,0.22],[0.4,0.43],[0.5,0.27],[0.43,0.34],[0.7,0.22],[0.62,0.22],[0.77,0.31],[0.77,0.4],[0.78,0.54],[0.77, 0.48],[0.25,0.335],[0.25,0.6],[0.25,0.51],[0.25,0.7],[0.25,0.78],[0.76,0.64],[0.777,0.73],[0.68,0.9],[0.755,0.86],[0.6,0.9],[0.42,0.9],[0.7,0.48],[0.59,0.55],[0.35,0.9],[0.25,0.9],[0.25,0.23],[0.49,0.9],[0.42,0.23],[0.32,0.63],[0.49,0.45],[0.25,0.42],[0.48,0.69]],
              [[0.34,0.22],[0.4,0.43],[0.5,0.27],[0.43,0.34],[0.7,0.22],[0.62,0.22],[0.77,0.31],[0.77,0.4],[0.78,0.54],[0.77, 0.48],[0.25,0.335],[0.25,0.6],[0.25,0.51],[0.25,0.7],[0.25,0.78],[0.76,0.64],[0.777,0.73],[0.68,0.9],[0.755,0.86],[0.6,0.9],[0.42,0.9],[0.7,0.48],[0.59,0.55],[0.35,0.9],[0.25,0.9],[0.25,0.23],[0.49,0.9],[0.42,0.23],[0.32,0.63],[0.49,0.45],[0.25,0.42],[0.48,0.69]],
              [[0.65,0.25],[0.48,0.4],[0.3,0.35],[0.25,0.55],[0.4,0.9],[0.25,0.84],[0.63,0.9],[0.7,0.65]],
              [[0.6,0.2],[0.45,0.2],[0.68,0.3],[0.75,0.3],[0.75,0.55],[0.75,0.7],[0.32,0.4],[0.38,0.42],[0.35,0.63],[0.43,0.63],[0.45,0.7],[0.5,0.7],[0.55,0.8],[0.45,0.81],[0.26,0.31],[0.3,0.23],[0.4,0.2],[0.7,0.2],[0.27,0.4],[0.29,0.6],[0.32,0.75],[0.32,0.9],[0.25,0.9],[0.35,0.82],[0.47,0.9],[0.63,0.8],[0.72,0.8],[0.72,0.9]],
              [[0.25, 0.9],[0.725,0.85],[0.25,0.5],[0.75,0.25],[0.75,0.65],[0.75,0.50],[0.4,0.27],[0.55,0.85],[0.53,0.27],[0.65,0.25],[0.25,0.3],[0.25,0.7]],
              [[0.6,0.4],[0.75,0.7],[0.5,0.72],[0.5,0.4],[0.6,0.55],[0.25,0.5],[0.4,0.27],[0.3,0.3],[0.25,0.66],[0.25, 0.9],[0.45,0.9],[0.65,0.9]],
              [[0.25,0.25],[0.25,0.48],[0.34,0.7],[0.51,0.45],[0.4,0.25],[0.7,0.25],[0.6,0.45],[0.53,0.7],[0.65,0.7],[0.75, 0.9]],
    )
    all_locs_all = (
                [[1,1],[4,4],[2,2],[1,1],[2,2],[2,2],[2,2],[2,2],[2,2],[2,2],[1,1],[4,4],[1,1],[4,4],[4,4],[2,2],[2,2],[3,3],[3,3],[3,3],[3,3],[2,2],[2,2],[3,3],[4,4],[1,1],[3,3],[1,1],[1,1],[2,2],[1,1],[4,4]],
                [[1,1],[4,4],[2,2],[1,1],[2,2],[2,2],[2,2],[2,2],[2,2],[2,2],[1,1],[4,4],[1,1],[4,4],[4,4],[2,2],[2,2],[3,3],[3,3],[3,3],[3,3]],
                [[1,4],[2,4],[2,3],[3,4],[2,3],[1,2],[1,4],[1,2]],
                [[2,3],[2,4],[1,2],[1,2],[1,3],[2,3],[1,2],[1,3],[1,3],[1,3],[1,3],[1,3],[1,3],[1,4],[2,3],[2,3],[2,3],[2,3],[2,3],[1,2],[1,2],[2,4],[3,4],[2,3],[2,3],[1,3],[1,3],[1,3]],
                [[1,3],[2,3],[2,4],[1,3],[2,3],[1,3],[1,2],[1,4],[1,2],[1,2],[2,4],[1,4]],
                [[1,3],[2,3],[2,3],[1,3],[2,3],[1,4],[1,2],[2,4],[1,3],[3,4],[2,4],[2,3]],
                [[2,4],[1,3],[1,3],[3,4],[2,4],[2,3],[1,4],[1,3],[1,3],[3,4]],
    )
    all_filenames = ('B3_insetB7.pdf' , 'B3_insetB6.pdf', 'B3_inset4.pdf', 'B3_inset5.pdf', 'B3_inset1.pdf', 'B3_inset2.pdf', 'B3_inset3.pdf')
    all_vmin = (-0.00001, -0.00001, -0.00001, -0.00001, -0.0001, -0.00001, -0.00001)
    all_fov_size = ((44,60,30,56,56,56,56)*u.arcsec).to(u.deg).value

    for sources, bboxes, locs_all, vmin, fov_size, filename in zip(all_sources, all_bboxes, all_locs_all, all_vmin, all_fov_size, all_filenames):
        print(filename)
        names = np.repeat(None, len(sources))

        zoomregions_auto = {}

        for i in range(len(sources)):
            key, reg_dict = create_zoomregion(sources[i], table, bboxes[i], locs=locs_all[i], name=names[i], vmin=vmin)
            zoomregions_auto[key] = reg_dict

        srcI_ind = np.where(table['D_ID'] == 30)[0]
        srcI_coord = coordinates.SkyCoord(table['RA_B3'][srcI_ind], table['DEC_B3'][srcI_ind], unit=u.degree)
        BL, TR = calc_bbox(srcI_coord, sidelength=fov_size)

        zoomregions_auto['(64,33)'] = {'bottomleft': coordinates.SkyCoord(ra=["5:35:14.435"],
                                                           dec=["-5:22:28.55"],
                                                           unit=(u.h, u.deg),
                                                           frame='icrs'),
                        'topright': coordinates.SkyCoord(ra=["5:35:14.403"],
                                                         dec=["-5:22:28.28"],
                                                         unit=(u.h, u.deg),
                                                         frame='icrs'),
                        'bbox':[0.57,0.66],'loc': 2,'l1':3,'l2':3,'min': -0.0001,'max': 0.0005,'zoom': 10}

        zoomregions_auto['(18,63)'] = {'topright': coordinates.SkyCoord(ra=["5:35:14.406"],
                                                           dec=["-5:22:33.162"],
                                                           unit=(u.h, u.deg),
                                                           frame='icrs'),
                        'bottomleft': coordinates.SkyCoord(ra=["5:35:14.445"],
                                                         dec=["-5:22:33.722"],
                                                         unit=(u.h, u.deg),
                                                         frame='icrs'),
                        'bbox':[0.6,0.4],'loc': 2,'l1':2,'l2':2,'min': -0.0001,'max': 0.0005,'zoom': 10}



        zoomregions = zoomregions_auto

        fn = 'Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'

        figure = inset_overlays(fn, zoomregions=zoomregions, directory=root,
                                vmin=-0.0005, vmax=0.003, bottomleft=BL, topright=TR,
                                make_overview=True,
                                savefn=filename.split(".")[0])

        figure = inset_overlays(fn, zoomregions=zoomregions, directory=root,
                                vmin=-0.0005, vmax=0.003, bottomleft=BL, topright=TR,
                                make_overview=False, colorbar=False,
                                savefn=filename.split(".")[0])

        figure = inset_overlays(fn, zoomregions=zoomregions, directory=root,
                                vmin=-0.0005, vmax=0.003, bottomleft=BL, topright=TR,
                                savefn=filename.split(".")[0])
