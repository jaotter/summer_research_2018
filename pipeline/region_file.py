from astrodendro import Dendrogram, pp_catalog
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import numpy as np
import radio_beam
import matplotlib.pyplot as pl

def compute_regions(min_val, min_del, npix, filename, reg_file, pdf=False):
    contfile = fits.open(filename)
    mywcs = wcs.WCS(contfile[0].header).celestial
    array = contfile[0].data.squeeze()
    beam = radio_beam.Beam.from_fits_header(contfile[0].header)
    d = Dendrogram.compute(array, min_value=min_val, min_delta = min_del, min_npix=npix, wcs=mywcs, verbose=False)
    
    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam
    pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    metadata['spatial_scale'] =  pixel_scale
    metadata['beam_major'] = beam.major
    metadata['beam_minor'] = beam.minor
    #metadata['wavelength'] = contfile[0].header['CRVAL3']*u.GHz
    metadata['wcs'] = mywcs
    cat = pp_catalog(d, metadata)
    
    with open(reg_file, 'w') as fh:
        fh.write("fk5\n")
        for row in cat:
            fh.write("ellipse({x_cen}, {y_cen}, {major_sigma}, "
                     "{minor_sigma}, {position_angle}) # text={{{_idx}}}\n"
    .format(**dict(zip(row.colnames, row))))
        
    if pdf == True:
        ax = pl.gca()
        ax.cla()
        pl.imshow(array, cmap='gray_r', interpolation='none', origin='lower',
            vmax=0.01, vmin=-0.001)
        pltr = d.plotter()
        for struct in d.leaves:
            cntr = pltr.plot_contour(ax, structure=struct, colors=['r'],
                                    linewidths=[0.9], zorder=5)
            if struct.parent:
                while struct.parent:
                    struct = struct.parent
                cntr_g = pltr.plot_contour(ax, structure=struct,
                                          colors=[(0,1,0,1)],
                                          linewidths=[0.2])
        pl.savefig(reg_file+'.pdf')
    return d, cat
