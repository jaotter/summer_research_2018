import numpy as np
import warnings
from astropy import units as u
from astropy.modeling import models, fitting
from radio_beam import Beam
import regions
from astropy import wcs
from astropy import coordinates
from astropy import log
from astropy.io import fits
from astropy.stats import mad_std
from astropy.utils.console import ProgressBar
import pylab as pl
import os
import signal, sys, time

def signal_handler(signal, frame):
    # this isn't strictly necessary, but I found that the loop wasn't
    # respecting my ctrl-C's, and I DEMAND my ctrl-C's be heard!!
    sys.exit(0)

STDDEV_TO_FWHM = np.sqrt(8*np.log(2))

def bbox_contains_bbox(bbox1,bbox2):
    """returns true if bbox2 is inside bbox1"""
    return ((bbox1.ixmax>bbox2.ixmax) & (bbox1.ixmin<bbox2.ixmin) &
            (bbox1.iymax>bbox2.iymax) & (bbox1.iymin<bbox2.iymin))

def sub_bbox_slice(bbox1, bbox2):
    """returns a slice from within bbox1 of bbox2"""
    if not bbox_contains_bbox(bbox1, bbox2):
        raise ValueError("bbox2 is not within bbox1")
    x0, dx = bbox2.ixmin-bbox1.ixmin, bbox2.ixmax-bbox2.ixmin
    y0, dy = bbox2.iymin-bbox1.iymin, bbox2.iymax-bbox2.iymin
    return (slice(y0, y0+dy), slice(x0, x0+dx),)

def slice_bbox_from_bbox(bbox1, bbox2):
    """
    Utility tool. Given two bboxes in the same coordinates, give the views of
    each box corresponding to the other.  For example, if you have an image
    ``im`` and two overlapping cutouts from that image ``cutout1`` and
    ``cutout2`` with bounding boxes ``bbox1`` and ``bbox2``, the returned views
    from this function give the regions ``cutout1[view1] = cutout2[view2]``
    """

    if bbox1.ixmin < bbox2.ixmin:
        blcx = bbox2.ixmin
    else:
        blcx = bbox1.ixmin
    if bbox1.ixmax > bbox2.ixmax:
        trcx = bbox2.ixmax
    else:
        trcx = bbox1.ixmax
    if bbox1.iymin < bbox2.iymin:
        blcy = bbox2.iymin
    else:
        blcy = bbox1.iymin
    if bbox1.iymax > bbox2.iymax:
        trcy = bbox2.iymax
    else:
        trcy = bbox1.iymax

    y0_1 = max(blcy-bbox1.iymin,0)
    x0_1 = max(blcx-bbox1.ixmin,0)
    y0_2 = max(blcy-bbox2.iymin,0)
    x0_2 = max(blcx-bbox2.ixmin,0)

    dy_1 = min(bbox1.iymax-blcy,trcy-blcy)
    dx_1 = min(bbox1.ixmax-blcx,trcx-blcx)
    dy_2 = min(bbox2.iymax-blcy,trcy-blcy)
    dx_2 = min(bbox2.ixmax-blcx,trcx-blcx)

    view1 = (slice(y0_1, y0_1+dy_1),
             slice(x0_1, x0_1+dx_1),)
    view2 = (slice(y0_2, y0_2+dy_2),
             slice(x0_2, x0_2+dx_2),)
    for slc in view1+view2:
        assert slc.start >= 0
        assert slc.stop >= 0
    return view1,view2




def bg_gaussfit(fitsfile, region, region_list,
                     radius=1.0*u.arcsec,
                     max_radius_in_beams=2,
                     max_offset_in_beams=1,
                     max_offset_in_beams_bg=10,
                     bg_stddev_x=40,
                     bg_stddev_y=40,
                     bg_mean_x=0,
                     bg_mean_y=0,
                     mask_size=1.5,
                     background_estimator=np.nanmedian,
                     noise_estimator=lambda x: mad_std(x, ignore_nan=True),
                     savepath=None,
                     prefix="",
                     covariance='param_cov',
                     raise_for_failure=False,
                    ):
    """
    Given a FITS filename and a region, fit a gaussian to the region
    with an input guess based on the beam size, and fit a gaussian to the background.

    Parameters
    ----------
    fitsfile : str
        Name of the FITS file
    region : region object
        single region from regions (see https://github.com/astropy/regions/)
    region_list : list of regions
        list of all regions - for masking out nearby sources
    radius : angular size
        The radius of the region around the region center to extract and
        include in the fit
    max_radius_in_beams : float
        The maximum allowed source radius in units of beam major axis
        (this is a limit passed to the fitter)
    max_offset_in_beams : float
        The maximum allowed offset of the source center from the guessed
        position
    max_offset_in_beams_bg : float
        same as above but for background gaussian
    bg_stddev_x : float
        Guess for standard deviation in x direction for background gaussian
    bg_stddev_y : float
        Same as above, in y direction
    bg_mean_x/y : float
        pixels away for background gaussian mean guess, with origin at center
    mask_size : float
        size in beams of mask to be applied to nearby sources
    background_estimator : function
        A function to apply to the background pixels (those not within 1 beam
        HWHM of the center) to estimate the background level.  The background
        will be subtracted before fitting.
    noise_estimator : function
        Function to apply to the whole data set to determine the noise level
        and therefore the appropriate per-pixel weight to get the correct
        normalization for the covariance matrix.
    savepath : str or None
        If specified, plots will be made and saved to this directory using the
        source name from the region metadata
    prefix : str
        The prefix to append to saved source names
    covariance : 'param_cov' or 'cov_x'
        Which covariance matrix should be used to estimate the parameter
        errors?  ``param_cov`` uses the diagonal of the reduced-chi^2-scaled
        covariance matrix to compute the parameter errors, while ``cov_x`` uses
        the unscaled errors.  See http://arxiv.org/abs/1009.2755 for a
        description, and criticism, of using the scaled covariance.
    raise_for_failure : bool
        If the fit was not successful, raise an exception
    """

    # need central coordinates of each object
    coords = coordinates.SkyCoord([reg.center for reg in region_list])

    fh = fits.open(fitsfile)
    data = fh[0].data.squeeze()
    header = fh[0].header
    datawcs = wcs.WCS(header).celestial
    beam = Beam.from_fits_header(header)
    pixscale = wcs.utils.proj_plane_pixel_area(datawcs)**0.5 * u.deg
    bmmin_px = (beam.minor.to(u.deg) / pixscale).decompose()
    bmmaj_px = (beam.major.to(u.deg) / pixscale).decompose()

    noise = noise_estimator(data)

    log.info("Noise estimate is {0} for file {1}".format(noise, fitsfile))

    fit_data = {}

    phot_reg = regions.CircleSkyRegion(center=region.center, radius=radius)
    pixreg = phot_reg.to_pixel(datawcs)
    mask = pixreg.to_mask()
    mask_cutout = mask.cutout(data)
    if mask_cutout is None:
        log.warning("The region failed to produce a cutout."
                        .format(reg))
        return null
    cutout = mask_cutout * mask.data
    cutout_mask = mask.data.astype('bool')

    smaller_phot_reg = regions.CircleSkyRegion(center=region.center,
                                                   radius=beam.major/2.) #FWHM->HWHM
    smaller_pixreg = smaller_phot_reg.to_pixel(datawcs)
    smaller_mask = smaller_pixreg.to_mask()
    smaller_cutout = smaller_mask.cutout(data) * smaller_mask.data

    # mask out (as zeros) neighboring sources within the fitting area
    srcind = None
    for ii, reg in enumerate(region_list):
        if reg.center.ra == region.center.ra and reg.center.dec == region.center.dec:
            srcind = ii
            
    nearby_matches = phot_reg.contains(coords, datawcs)
    if any(nearby_matches):
        inds = np.where(nearby_matches)[0].tolist()
        inds.remove(srcind)
        for ind in inds:
            maskoutreg = regions.EllipseSkyRegion(center=region_list[ind].center,
                                                      width=mask_size*beam.major,
                                                      height=mask_size*beam.minor,
                                                      angle=beam.pa+90*u.deg,
                                                     )
            mpixreg = maskoutreg.to_pixel(datawcs)
            mmask = mpixreg.to_mask()

            view, mview = slice_bbox_from_bbox(mask.bbox, mmask.bbox)
            cutout_mask[view] &= ~mmask.data.astype('bool')[mview]
            cutout = cutout * cutout_mask


    background_mask = cutout_mask.copy().astype('bool')
    background_mask[sub_bbox_slice(mask.bbox, smaller_mask.bbox)] &= ~smaller_mask.data.astype('bool')
    background = background_estimator(cutout[background_mask])

    sz = cutout.shape[0]
    mx = np.nanmax(smaller_cutout)
    ampguess = mx-background
    imtofit = np.nan_to_num((cutout-background)*mask.data)
    src_gauss = [ampguess, sz/2, bmmaj_px.value, bmmin_px.value, beam.pa.value]
    bg_gauss = [background, bg_mean_x + sz/2, bg_mean_y + sz/2, bg_stddev_x, bg_stddev_y, beam.pa.value]
    bnds = [max_radius_in_beams, max_offset_in_beams, max_offset_in_beams_bg]

    result, fit_info, chi2, fitter = gaussfit_image(image=imtofit,
                                                    gauss_params=src_gauss,
                                                    bg_gauss_params=bg_gauss,
                                                    bound_params=bnds,
                                                    weights=1/noise**2,
                                                    plot=savepath is not None,
    )
    sourcename = region.meta['text'].strip('{}')

    if savepath is not None:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            bmarr = beam.as_kernel(pixscale=pixscale, x_size=sz, y_size=sz).array
            assert bmarr.max() > 0
            bm_ellipse = beam.ellipse_to_plot(sz/2, sz/2., pixscale)
            bm_ellipse.set_facecolor('none')
            bm_ellipse.set_edgecolor('r')
            pl.gca().add_patch(bm_ellipse)
            #pl.contour(bmarr, levels=[0.317*bmarr.max()], colors=['r'])
        pl.savefig(os.path.join(savepath, '{0}{1}.png'.format(prefix, sourcename)),
                   bbox_inches='tight')

    if covariance not in fit_info or fit_info[covariance] is None:
        fit_info[covariance] = np.zeros([6,6])
        success = False
    else:
        success = True

    cx,cy = pixreg.bounding_box.ixmin+result.x_mean_0, pixreg.bounding_box.iymin+result.y_mean_0
    clon,clat = datawcs.wcs_pix2world(cx, cy, 0)

    major,minor = (result.x_stddev_0 * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                   result.y_stddev_0 * STDDEV_TO_FWHM * pixscale.to(u.arcsec))
    majind, minind = 3,4
    pa = (result.theta_0*u.rad).to(u.deg)
    if minor > major:
        major,minor = minor,major
        majind,minind = minind,majind
        pa += 90*u.deg

    fitted_gaussian_as_beam = Beam(major=major, minor=minor, pa=pa)
    try:
        deconv_fit = fitted_gaussian_as_beam.deconvolve(beam)
        deconv_major, deconv_minor, deconv_pa = (deconv_fit.major,
                                                 deconv_fit.minor,
                                                 deconv_fit.pa)
        deconv_maj_err = fit_info[covariance][majind,majind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec)
        deconv_min_err = fit_info[covariance][minind,minind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec)
        deconv_pa_err = fit_info[covariance][5,5]**0.5 * u.deg
    except ValueError:
        print("Could not deconvolve {0} from {1}".format(beam.__repr__(), fitted_gaussian_as_beam.__repr__()))
        deconv_major, deconv_minor, deconv_pa = np.nan, np.nan, np.nan
        deconv_maj_err, deconv_min_err, deconv_pa_err = np.nan, np.nan, np.nan
        fit_data[sourcename] = {'amplitude': result.amplitude_0,
                                'center_x': float(clon)*u.deg,
                                'center_y': float(clat)*u.deg,
                                'fwhm_major': major,
                                'fwhm_minor': minor,
                                'pa': pa,
                                'deconv_fwhm_major': deconv_major,
                                'e_deconv_fwhm_major' : deconv_maj_err,
                                'deconv_fwhm_minor': deconv_minor,
                                'e_deconv_fwhm_minor': deconv_min_err,
                                'deconv_pa': deconv_pa,
                                'e_deconv_pa': deconv_pa_err,
                                'chi2': chi2,
                                'chi2/n': chi2/mask.data.sum(),
                                'e_amplitude': fit_info[covariance][0,0]**0.5,
                                'e_center_x': fit_info[covariance][1,1]**0.5*pixscale,
                                'e_center_y': fit_info[covariance][2,2]**0.5*pixscale,
                                'e_fwhm_major': fit_info[covariance][majind,majind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                                'e_fwhm_minor': fit_info[covariance][minind,minind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                                'e_pa': fit_info[covariance][5,5]**0.5 * u.deg,
                                'success': success,
                                'ampguess': ampguess,
                                'peak': mx,
                                'fit_info': fit_info,
        }

    if raise_for_failure and not success:
        raise ValueError("Fit failed.")

    signal.signal(signal.SIGINT, signal_handler)

    return fit_data


def gaussfit_image(image, gauss_params, bg_gauss_params, bound_params, weights=None,
                   fitter=fitting.LevMarLSQFitter(), plot=False):
    """
    Fit a gaussian to an image and optionally plot the data, the fitted
    gaussian, and the residual.

    Parameters
    ----------
    image : 2-dimensional array
        The image to be fit.  Cannot contain any NaNs.
    gaussian : list of floats
        [amplitude guess, x and y mean (center), beam major pix, beam minor pix, beam pa]
    bg_gaussian : list of floats
        [background guess, bg mean x , bg mean y, bg_stddev_x, bg_stddev_y, beam pa]
    fitter : `astropy.modeling.fitting.Fitter`
        A fitter instance.  Can be any of the optimizers, in principle, but it
        needs to take keywords ``weight`` and ``maxiter``.
    plot : bool
        Make the "diagnostic plot" showing the image, the fitted image, the
        residual, and the image with the fits contoured over it?

    Returns
    -------
    fitted : `astropy.modeling.fitting.Fitter`
        The fitter instance
    fitter.fit_info : dict
        The dictionary containing the ``fit_info`` from the fitter
    residualsquaredsum : float
        The sum of the squares of the residual, e.g., chi^2.
    """
    yy, xx = np.mgrid[:image.shape[0], :image.shape[1]]

    #print(gauss_params)
    #print(bg_gauss_params)
    src_gaussian = models.Gaussian2D(amplitude=gauss_params[0],
                                   x_mean=gauss_params[1],
                                   y_mean=gauss_params[1],
                                   x_stddev=gauss_params[2]/STDDEV_TO_FWHM,
                                   y_stddev=gauss_params[3]/STDDEV_TO_FWHM,
                                   theta=gauss_params[4],
                                   bounds={'x_stddev':(gauss_params[3]/STDDEV_TO_FWHM*0.75,
                                                       gauss_params[2]*bound_params[0]/STDDEV_TO_FWHM),
                                           'y_stddev':(gauss_params[3]/STDDEV_TO_FWHM*0.75,
                                                       gauss_params[2]*bound_params[0]/STDDEV_TO_FWHM),
                                           'x_mean':(gauss_params[1]-bound_params[1]*gauss_params[2]/STDDEV_TO_FWHM,
                                                     gauss_params[1]+bound_params[1]*gauss_params[2]/STDDEV_TO_FWHM),
                                           'y_mean':(gauss_params[1]-bound_params[1]*gauss_params[2]/STDDEV_TO_FWHM,
                                                     gauss_params[1]+bound_params[1]*gauss_params[2]/STDDEV_TO_FWHM),
                                           'amplitude':(gauss_params[0]*0.9, gauss_params[0]*1.1)
                                          }
                                     )
    bg_gaussian = models.Gaussian2D(amplitude=0.1*gauss_params[0],
                                   x_mean=bg_gauss_params[1],
                                   y_mean=bg_gauss_params[2],
                                   x_stddev=bg_gauss_params[3]/STDDEV_TO_FWHM,
                                   y_stddev=bg_gauss_params[4]/STDDEV_TO_FWHM,
                                   theta=bg_gauss_params[5],
                                   bounds={'amplitude':(bg_gauss_params[0]*0.01, 0.5*gauss_params[0]),
                                           'x_mean':(bg_gauss_params[1]-bound_params[2]*gauss_params[3]/STDDEV_TO_FWHM,
                                                     bg_gauss_params[1]+bound_params[2]*gauss_params[3]/STDDEV_TO_FWHM),
                                           'y_mean':(bg_gauss_params[2]-bound_params[2]*gauss_params[3]/STDDEV_TO_FWHM,
                                                     bg_gauss_params[2]+bound_params[2]*gauss_params[3]/STDDEV_TO_FWHM)}
    )
    
    gauss_init = src_gaussian + bg_gaussian

    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        fitted = fitter(gauss_init, xx, yy, image, weights=weights,
                        maxiter=1000)

    print(fitted)
    fitim = fitted(xx,yy)
    fitim_src = fitted[0](xx,yy)
    residual = image-fitim
    residualsquaredsum = np.nansum(residual**2*weights)

    if plot:
        pl.clf()
        ax1 = pl.subplot(2,3,1)
        im = ax1.imshow(image, cmap='viridis', origin='lower',
                        interpolation='nearest')
        vmin, vmax = im.get_clim()
        ax2 = pl.subplot(2,3,2)
        ax2.imshow(fitim, cmap='viridis', origin='lower',
                   interpolation='nearest', vmin=vmin, vmax=vmax)
        ax3 = pl.subplot(2,3,3)
        ax3.imshow(fitim_src, cmap='viridis', origin='lower',
                   interpolation='nearest', vmin=vmin, vmax=vmax)
        ax4 = pl.subplot(2,3,4)
        ax4.imshow(residual, cmap='viridis', origin='lower',
                   interpolation='nearest', vmin=vmin, vmax=vmax)
        ax5 = pl.subplot(2,3,5)
        im = ax5.imshow(image, cmap='viridis', origin='lower',
                        interpolation='nearest')
        ax6 = pl.subplot(2,3,6)
        
        vmin, vmax = im.get_clim()
        scalefactor = fitim.max()
        if scalefactor < 0:
            scalefactor = fitim.max() - fitim.min()
        ax5.contour(fitim, levels=np.array([0.00269, 0.0455, 0.317])*scalefactor,
                    colors=['w']*4)
        axlims = ax4.axis()
        ax5.plot(fitted.x_mean_0, fitted.y_mean_0, 'w+')
        ax5.axis(axlims)

        ax6 = pl.subplot(2,3,6)
        im = ax6.imshow(image, cmap='viridis', origin='lower',
                        interpolation='nearest')
        ax6.contour(fitim_src, levels=np.array([0.00269, 0.0455, 0.317])*scalefactor,
                    colors=['w']*4)
        axlims = ax6.axis()
        ax6.plot(fitted.x_mean_0, fitted.y_mean_0, 'w+')
        ax6.axis(axlims)
        
    return fitted, fitter.fit_info, residualsquaredsum, fitter


def gaussfit_sub_bg(fitsfile, region, region_list,
                     bg_stddev_x,
                     bg_stddev_y,
                     bg_mean_x,
                     bg_mean_y,
                     bg_theta,
                     bg_amp,
                     radius=1.0*u.arcsec,
                     max_radius_in_beams=2,
                     max_offset_in_beams=1,
                     mask_size=1.5,
                     background_estimator=np.nanmedian,
                     noise_estimator=lambda x: mad_std(x, ignore_nan=True),
                     savepath=None,
                     prefix="",
                     covariance='param_cov',
                     raise_for_failure=False,
                    ):
    """
    Given a FITS filename and a region, fit a gaussian to the region
    with an input guess based on the beam size, and fit a gaussian to the background.

    Parameters
    ----------
    fitsfile : str
        Name of the FITS file
    region : region object
        single region from regions (see https://github.com/astropy/regions/)
    region_list : list of regions
        list of all regions - for masking out nearby sources
    radius : angular size
        The radius of the region around the region center to extract and
        include in the fit
    max_radius_in_beams : float
        The maximum allowed source radius in units of beam major axis
        (this is a limit passed to the fitter)
    max_offset_in_beams : float
        The maximum allowed offset of the source center from the guessed
        position
    max_offset_in_beams_bg : float
        same as above but for background gaussian
    bg_stddev_x : float
        Standard deviation in x direction for background gaussian
    bg_stddev_y : float
        Same as above, in y direction
    bg_mean_x/y : float
        pixels away for background gaussian mean, with origin at center
    bg_theta : float
        position angle of background gaussian
    bg_amp : float
        background amplitude
    mask_size : float
        size in beams of mask to be applied to nearby sources
    background_estimator : function
        A function to apply to the background pixels (those not within 1 beam
        HWHM of the center) to estimate the background level.  The background
        will be subtracted before fitting.
    noise_estimator : function
        Function to apply to the whole data set to determine the noise level
        and therefore the appropriate per-pixel weight to get the correct
        normalization for the covariance matrix.
    savepath : str or None
        If specified, plots will be made and saved to this directory using the
        source name from the region metadata
    prefix : str
        The prefix to append to saved source names
    covariance : 'param_cov' or 'cov_x'
        Which covariance matrix should be used to estimate the parameter
        errors?  ``param_cov`` uses the diagonal of the reduced-chi^2-scaled
        covariance matrix to compute the parameter errors, while ``cov_x`` uses
        the unscaled errors.  See http://arxiv.org/abs/1009.2755 for a
        description, and criticism, of using the scaled covariance.
    raise_for_failure : bool
        If the fit was not successful, raise an exception
    """

    # need central coordinates of each object
    coords = coordinates.SkyCoord([reg.center for reg in region_list])

    fh = fits.open(fitsfile)
    data = fh[0].data.squeeze()
    header = fh[0].header
    datawcs = wcs.WCS(header).celestial
    beam = Beam.from_fits_header(header)
    pixscale = wcs.utils.proj_plane_pixel_area(datawcs)**0.5 * u.deg
    bmmin_px = (beam.minor.to(u.deg) / pixscale).decompose()
    bmmaj_px = (beam.major.to(u.deg) / pixscale).decompose()

    noise = noise_estimator(data)

    log.info("Noise estimate is {0} for file {1}".format(noise, fitsfile))

    fit_data = {}

    phot_reg = regions.CircleSkyRegion(center=region.center, radius=radius)
    pixreg = phot_reg.to_pixel(datawcs)
    mask = pixreg.to_mask()
    mask_cutout = mask.cutout(data)
    if mask_cutout is None:
        log.warning("The region failed to produce a cutout."
                        .format(reg))
        return null
    cutout = mask_cutout * mask.data
    cutout_mask = mask.data.astype('bool')

    smaller_phot_reg = regions.CircleSkyRegion(center=region.center,
                                                   radius=beam.major/2.) #FWHM->HWHM
    smaller_pixreg = smaller_phot_reg.to_pixel(datawcs)
    smaller_mask = smaller_pixreg.to_mask()
    smaller_cutout = smaller_mask.cutout(data) * smaller_mask.data

    # mask out (as zeros) neighboring sources within the fitting area
    srcind = None
    for ii, reg in enumerate(region_list):
        if reg.center.ra == region.center.ra and reg.center.dec == region.center.dec:
            srcind = ii
    print(srcind)
    nearby_matches = phot_reg.contains(coords, datawcs)
    if any(nearby_matches):
        inds = np.where(nearby_matches)[0].tolist()
        inds.remove(srcind)
        for ind in inds:
            maskoutreg = regions.EllipseSkyRegion(center=region_list[ind].center,
                                                      width=mask_size*beam.major,
                                                      height=mask_size*beam.minor,
                                                      angle=beam.pa+90*u.deg,
                                                     )
            mpixreg = maskoutreg.to_pixel(datawcs)
            mmask = mpixreg.to_mask()

            view, mview = slice_bbox_from_bbox(mask.bbox, mmask.bbox)
            cutout_mask[view] &= ~mmask.data.astype('bool')[mview]
            cutout = cutout * cutout_mask


        background_mask = cutout_mask.copy().astype('bool')
        background_mask[sub_bbox_slice(mask.bbox, smaller_mask.bbox)] &= ~smaller_mask.data.astype('bool')
        background = background_estimator(cutout[background_mask])

        sz = cutout.shape[0]
        mx = np.nanmax(smaller_cutout)
        ampguess = mx-background
        imtofit = np.nan_to_num((cutout-background)*mask.data)
        src_gauss = [ampguess, sz/2, bmmaj_px.value, bmmin_px.value, beam.pa.value]
        bg_gauss = [background, bg_mean_x + sz/2, bg_mean_y + sz/2, bg_stddev_x, bg_stddev_y, beam.pa.value]
        bnds = [max_radius_in_beams, max_offset_in_beams, max_offset_in_beams_bg]
        result, fit_info, chi2, fitter = gaussfit_image(image=imtofit,
                                                        gauss_params=src_gauss,
                                                        bg_gauss_params=bg_gauss,
                                                        bound_params=bnds,
                                                        weights=1/noise**2,
                                                        plot=savepath is not None,
                                                       )
        sourcename = region.meta['text'].strip('{}')

        if savepath is not None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                bmarr = beam.as_kernel(pixscale=pixscale, x_size=sz, y_size=sz).array
            assert bmarr.max() > 0
            bm_ellipse = beam.ellipse_to_plot(sz/2, sz/2., pixscale)
            bm_ellipse.set_facecolor('none')
            bm_ellipse.set_edgecolor('r')
            pl.gca().add_patch(bm_ellipse)
            #pl.contour(bmarr, levels=[0.317*bmarr.max()], colors=['r'])
            pl.savefig(os.path.join(savepath, '{0}{1}.png'.format(prefix, sourcename)),
                       bbox_inches='tight')

        if covariance not in fit_info or fit_info[covariance] is None:
            fit_info[covariance] = np.zeros([6,6])
            success = False
        else:
            success = True

        cx,cy = pixreg.bounding_box.ixmin+result.x_mean_0, pixreg.bounding_box.iymin+result.y_mean_0
        clon,clat = datawcs.wcs_pix2world(cx, cy, 0)

        major,minor = (result.x_stddev_0 * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                       result.y_stddev_0 * STDDEV_TO_FWHM * pixscale.to(u.arcsec))
        majind, minind = 3,4
        pa = (result.theta_0*u.rad).to(u.deg)
        if minor > major:
            major,minor = minor,major
            majind,minind = minind,majind
            pa += 90*u.deg

        fitted_gaussian_as_beam = Beam(major=major, minor=minor, pa=pa)
        try:
            deconv_fit = fitted_gaussian_as_beam.deconvolve(beam)
            deconv_major, deconv_minor, deconv_pa = (deconv_fit.major,
                                                     deconv_fit.minor,
                                                     deconv_fit.pa)
            deconv_maj_err = fit_info[covariance][majind,majind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec)
            deconv_min_err = fit_info[covariance][minind,minind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec)
            deconv_pa_err = fit_info[covariance][5,5]**0.5 * u.deg
        except ValueError:
            print("Could not deconvolve {0} from {1}".format(beam.__repr__(), fitted_gaussian_as_beam.__repr__()))
            deconv_major, deconv_minor, deconv_pa = np.nan, np.nan, np.nan
            deconv_maj_err, deconv_min_err, deconv_pa_err = np.nan, np.nan, np.nan
        fit_data[sourcename] = {'amplitude': result.amplitude_0,
                                'center_x': float(clon)*u.deg,
                                'center_y': float(clat)*u.deg,
                                'fwhm_major': major,
                                'fwhm_minor': minor,
                                'pa': pa,
                                'deconv_fwhm_major': deconv_major,
                                'e_deconv_fwhm_major' : deconv_maj_err,
                                'deconv_fwhm_minor': deconv_minor,
                                'e_deconv_fwhm_minor': deconv_min_err,
                                'deconv_pa': deconv_pa,
                                'e_deconv_pa': deconv_pa_err,
                                'chi2': chi2,
                                'chi2/n': chi2/mask.data.sum(),
                                'e_amplitude': fit_info[covariance][0,0]**0.5,
                                'e_center_x': fit_info[covariance][1,1]**0.5*pixscale,
                                'e_center_y': fit_info[covariance][2,2]**0.5*pixscale,
                                'e_fwhm_major': fit_info[covariance][majind,majind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                                'e_fwhm_minor': fit_info[covariance][minind,minind]**0.5 * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                                'e_pa': fit_info[covariance][5,5]**0.5 * u.deg,
                                'success': success,
                                'ampguess': ampguess,
                                'peak': mx,
                                'fit_info': fit_info,
                               }

        if raise_for_failure and not success:
            raise ValueError("Fit failed.")

        signal.signal(signal.SIGINT, signal_handler)

    return fit_data
