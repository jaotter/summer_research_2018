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




def gaussfit_catalog(fitsfile, region_list, radius=1.0*u.arcsec,
                     max_radius_in_beams=2,
                     max_offset_in_beams=1,
                     background_estimator=np.nanmedian,
                     noise_estimator=lambda x: mad_std(x, ignore_nan=True),
                     savepath=None,
                     prefix="",
                     covariance='param_cov',
                     raise_for_failure=False,
                    ):
    """
    Given a FITS filename and a list of regions, fit a gaussian to each region
    with an input guess based on the beam size.

    Parameters
    ----------
    fitsfile : str
        Name of the FITS file
    region_list : list
        List of regions (see https://github.com/astropy/regions/)
    radius : angular size
        The radius of the region around the region center to extract and
        include in the fit
    max_radius_in_beams : float
        The maximum allowed source radius in units of beam major axis
        (this is a limit passed to the fitter)
    max_offset_in_beams : float
        The maximum allowed offset of the source center from the guessed
        position
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

    pb = ProgressBar(len(region_list))

    for ii,reg in enumerate(region_list):

        phot_reg = regions.CircleSkyRegion(center=reg.center, radius=radius)
        pixreg = phot_reg.to_pixel(datawcs)
        mask = pixreg.to_mask()
        mask_cutout = mask.cutout(data)
        if mask_cutout is None:
            log.warning("Skipping region {0} because it failed to produce a cutout."
                        .format(reg))
            continue
        cutout = mask_cutout * mask.data
        cutout_mask = mask.data.astype('bool')

        smaller_phot_reg = regions.CircleSkyRegion(center=reg.center,
                                                   radius=beam.major/2.) #FWHM->HWHM
        smaller_pixreg = smaller_phot_reg.to_pixel(datawcs)
        smaller_mask = smaller_pixreg.to_mask()
        smaller_cutout = smaller_mask.cutout(data) * smaller_mask.data

        # mask out (as zeros) neighboring sources within the fitting area
        nearby_matches = phot_reg.contains(coords, datawcs)
        if any(nearby_matches):
            inds = np.where(nearby_matches)[0].tolist()
            inds.remove(ii)
            for ind in inds:
                maskoutreg = regions.EllipseSkyRegion(center=region_list[ind].center,
                                                      width=beam.major,
                                                      height=beam.minor,
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

        p_init = models.Gaussian2D(amplitude=ampguess,
                                   x_mean=sz/2,
                                   y_mean=sz/2,
                                   x_stddev=bmmaj_px/STDDEV_TO_FWHM,
                                   y_stddev=bmmin_px/STDDEV_TO_FWHM,
                                   theta=beam.pa,
                                   bounds={'x_stddev':(bmmin_px/STDDEV_TO_FWHM*0.75,
                                                       bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                           'y_stddev':(bmmin_px/STDDEV_TO_FWHM*0.75,
                                                       bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                           'x_mean':(sz/2-max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM,
                                                     sz/2+max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM),
                                           'y_mean':(sz/2-max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM,
                                                     sz/2+max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM),
                                           'amplitude':(ampguess*0.9, ampguess*1.1)
                                          }
                                  )

        imtofit = np.nan_to_num((cutout-background)*mask.data)
        result, fit_info, chi2, fitter = gaussfit_image(image=imtofit,
                                                        gaussian=p_init,
                                                        weights=1/noise**2,
                                                        plot=savepath is not None,
                                                       )
        sourcename = reg.meta['text'].strip('{}')

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

        cx,cy = pixreg.bounding_box.ixmin+result.x_mean, pixreg.bounding_box.iymin+result.y_mean
        clon,clat = datawcs.wcs_pix2world(cx, cy, 0)

        major,minor = (result.x_stddev * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                       result.y_stddev * STDDEV_TO_FWHM * pixscale.to(u.arcsec))
        majind, minind = 3,4
        pa = (result.theta*u.rad).to(u.deg)
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
        except ValueError:
            print("Could not deconvolve {0} from {1}".format(beam.__repr__(), fitted_gaussian_as_beam.__repr__()))
            deconv_major, deconv_minor, deconv_pa = np.nan, np.nan, np.nan

        fit_data[sourcename] = {'amplitude': result.amplitude,
                                'center_x': float(clon)*u.deg,
                                'center_y': float(clat)*u.deg,
                                'fwhm_major': major,
                                'fwhm_minor': minor,
                                'pa': pa,
                                'deconv_fwhm_major': deconv_major,
                                'deconv_fwhm_minor': deconv_minor,
                                'deconv_pa': deconv_pa,
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

        pb.update(ii)
        signal.signal(signal.SIGINT, signal_handler)

    return fit_data


def gaussfit_image(image, gaussian, weights=None,
                   fitter=fitting.LevMarLSQFitter(), plot=False):
    """
    Fit a gaussian to an image and optionally plot the data, the fitted
    gaussian, and the residual.

    Parameters
    ----------
    image : 2-dimensional array
        The image to be fit.  Cannot contain any NaNs.  Should have zero
        background, since the model (if it's a Gaussian model) does not include
        a background.
    gaussian : `astropy.modeling.Model`
        An astropy model object with guesses included.  Given the name of this
        function, it should be a `~astropy.models.Gaussian2D` model, but
        technically it can be any model.
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
    
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        fitted = fitter(gaussian, xx, yy, image, weights=weights,
                        maxiter=1000)

    fitim = fitted(xx,yy)
    residual = image-fitim
    residualsquaredsum = np.nansum(residual**2*weights)

    if plot:
        pl.clf()
        ax1 = pl.subplot(2,2,1)
        im = ax1.imshow(image, cmap='viridis', origin='lower',
                        interpolation='nearest')
        vmin, vmax = im.get_clim()
        ax2 = pl.subplot(2,2,2)
        ax2.imshow(fitim, cmap='viridis', origin='lower',
                   interpolation='nearest', vmin=vmin, vmax=vmax)
        ax3 = pl.subplot(2,2,3)
        ax3.imshow(residual, cmap='viridis', origin='lower',
                   interpolation='nearest', vmin=vmin, vmax=vmax)
        ax4 = pl.subplot(2,2,4)
        im = ax4.imshow(image, cmap='viridis', origin='lower',
                        interpolation='nearest')
        vmin, vmax = im.get_clim()
        scalefactor = fitim.max()
        if scalefactor < 0:
            scalefactor = fitim.max() - fitim.min()
        ax4.contour(fitim, levels=np.array([0.00269, 0.0455, 0.317])*scalefactor,
                    colors=['w']*4)
        axlims = ax4.axis()
        ax4.plot(fitted.x_mean, fitted.y_mean, 'w+')
        ax4.axis(axlims)
    
    return fitted, fitter.fit_info, residualsquaredsum, fitter
