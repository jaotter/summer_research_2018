from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.io import fits
from scipy.stats import median_abs_deviation
from latex_info import rounded
import radio_beam
import regions
import numpy as np
import astropy.units as u

def mask(reg, cutout):
    #masks everything except the region  
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')

def measure_rms(coord, data, img_wcs, annulus_radius=0.1*u.arcsecond):
    
    pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

    annulus_width = 15 #pix


    center_coord_pix = coord.to_pixel(img_wcs)
    cutout = Cutout2D(data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
    cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
        
    innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value)
    outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value+annulus_width)

    innerann_mask = innerann_reg.to_mask()

    annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)

    # Calculate the SNR and aperture flux sums                                                                                                                                    
    pixels_in_annulus = cutout.data[annulus_mask.astype('bool')]
    bg_rms = median_abs_deviation(pixels_in_annulus)

    return bg_rms

def nondet_table():
    IRtab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_may21_full_edit.fits')
    MLLA = Table.read('/home/jotter/nrao/tables/MLLA_02_IR.fit')

    b3_rad = 137.6*u.arcsecond / 2
    b3_center = SkyCoord(ra='05:35:14.51',dec='-05:22:30.56', unit=(u.hourangle, u.degree))
    
    MLLA_coord = SkyCoord(MLLA['RAJ2000'], MLLA['DEJ2000'], unit=u.degree)
    dist = MLLA_coord.separation(b3_center)

    fov_ind = np.where(dist < b3_rad)
    
    MLLA_fov = MLLA[fov_ind]

    print(len(MLLA_fov), len(IRtab))

    MLLA_det_ind = []
    for row in IRtab:
        mlla_ind = np.where(MLLA_fov['MLLA'] == row['MLLA'])[0]
        if len(mlla_ind) > 1:
            mlla_ind2 = np.where(MLLA_fov[mlla_ind]['m_MLLA'] == row['m_MLLA'])[0]
            mlla_ind = mlla_ind[mlla_ind2]#MLLA_fov[mlla_ind]
        MLLA_det_ind.append(mlla_ind[0])

    MLLA_fov_nondet_ind = np.setdiff1d(np.arange(len(MLLA_fov)), np.array(MLLA_det_ind))
    MLLA_nondet = MLLA_fov[MLLA_fov_nondet_ind]

    MLLA_nondet.write('/home/jotter/nrao/summer_research_2018/tables/IR_nondet_may21_full.fits', overwrite=True)


def mlla_nondet():
    MLLA_nondet = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_nondet_may21_full.fits')

    b3_fl = fits.open('/home/jotter/nrao/images/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.huge.deepmask.image.tt0.pbcor.fits')
    header = b3_fl[0].header
    img_wcs = WCS(header)
    data = b3_fl[0].data
    beam = radio_beam.Beam.from_fits_header(header)
    pixel_scale = np.abs(img_wcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    #ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
    
    MLLA_coord = SkyCoord(ra=MLLA_nondet['RAJ2000'], dec=MLLA_nondet['DEJ2000'], unit=u.degree)
    annulus_radius = 0.1*u.arcsecond
    annulus_radius_pix = (annulus_radius.to(u.degree)/pixel_scale).decompose()

    annulus_width = 15 #pix
    
    ulim_fluxes = []
    
    for ind in range(len(MLLA_coord)):
        center_coord = MLLA_coord[ind]
        center_coord_pix = center_coord.to_pixel(img_wcs)
        cutout = Cutout2D(data, center_coord_pix, annulus_radius*2.5, img_wcs, mode='partial')
        cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
        
        innerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value)
        outerann_reg = regions.CirclePixelRegion(cutout_center, annulus_radius_pix.value+annulus_width)

        innerann_mask = innerann_reg.to_mask()
        
        annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)
        
        # Calculate the SNR and aperture flux sums                                                                                                                                    
        pixels_in_annulus = cutout.data[annulus_mask.astype('bool')]
        bg_rms = median_abs_deviation(pixels_in_annulus)

        ulim_fluxes.append(3*bg_rms)

    ulim_fluxes = np.array(ulim_fluxes)*1000*u.mJy #in mJy
    MLLA_nondet['B3_flux_ulim'] = ulim_fluxes
    MLLA_nondet.write('/home/jotter/nrao/summer_research_2018/tables/IR_nondet_may21_full_ulim.fits', overwrite=True)
    

def b3_nondet(srcs_ID, band):
    if band == 'B6':
        img_path = '/home/jotter/nrao/images/Orion_SourceI_B6_continuum_r0.5.clean1mJy.150mplus.huge.image.tt0.pbcor.fits'
    if band == 'B7':
        img_path = '/home/jotter/nrao/images/Orion_SourceI_B7_continuum_r0.5.clean0.05mJy.250klplus.deepmask.image.tt0.pbcor.fits'
        
    tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits')

    fl = fits.open(img_path)
    header = fl[0].header
    img_wcs = WCS(header)
    data = fl[0].data
    beam = radio_beam.Beam.from_fits_header(header)
    
    coord_all = SkyCoord(ra=tab['RA_B3'][srcs_ID], dec=tab['DEC_B3'][srcs_ID], unit=u.degree)
    
    ulim_fluxes = []
    
    for coord in coord_all:
        bg_rms = measure_rms(coord, data, img_wcs)
        ulim_fluxes.append(3*bg_rms)
    print(ulim_fluxes)
    ulim_fluxes = np.array(ulim_fluxes)*1000*u.mJy #in mJy
    band_ulim_flux = np.repeat(np.nan, len(tab))
    band_ulim_flux[srcs_ID] = ulim_fluxes
    tab[f'{band}_flux_ulim'] = band_ulim_flux
    tab.write('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_may21_ulim_mask.fits', overwrite=True)


def paper_table():
    tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_nondet_may21_full_ulim_hc2000.fits')

    #MLLA = Column(np.array(np.repeat('-', len(data)), dtype='S10'), name='MLLA')
    MLLA_col = []
    b3_ulim = []
    altid = []
    for ir_row in tab:
        if ir_row['m_MLLA'] == ' ':
            mlla = str(ir_row['MLLA'])
        else:
            mlla = str(ir_row['MLLA']) + str(ir_row['m_MLLA']).lower()
        MLLA_col.append(mlla)
        b3_ulim_val, val = rounded(ir_row['B3_flux_ulim'], ir_row['B3_flux_ulim']/10, extra=0)
        b3_ulim.append(b3_ulim_val)
        if len(ir_row['AltID']) > 1:
            altid.append(ir_row['AltID'])
        else:
            altid.append('')
        
    paper_tab = Table((MLLA_col, altid,  b3_ulim), names=('MLLA', 'Alternate ID', '3mm Flux Upper Limit'), units=(None, None, u.mJy))
    
    print(paper_tab)
    paper_tab.write('/home/jotter/nrao/tables/latex_tables/final_tables/table2_full.fits', overwrite=True)
    paper_tab.write('/home/jotter/nrao/tables/latex_tables/table2_full.txt', format='latex', overwrite=True)
                      
#nondet_table()
#mlla_nondet()
#paper_table()
b6_nondet_srcs = [3, 14, 19, 37, 38, 43, 53, 56, 58]
b7_nondet_srcs = [37, 38]

b3_nondet(b6_nondet_srcs, 'B6')
b3_nondet(b7_nondet_srcs, 'B7')
