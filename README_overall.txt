/pipeline/ contains files relating to source detection and flux measurements

/flux_plots/ is for making SEDs and other flux plots

/proper_motions/ has scripts for calculating proper motions

As for the data tables, they are all located in
/lustre/aoc/students/jotter/dendro_catalogs/
The 'master' table: 'master_500klplus_B3_ref.fits' contains the fluxes in the r-2.clean0.1mJy.500klplus.deepmask images for B3, B6, and B7 using apertures from the original B3 gaussian fits. Note that the values of the FWHM major and minor axes, and the position angle are from gaussian fits from these images. Also note the positions are still from the original B3 data (which maybe should be changed). This table also has measurements from the ALMA 340GHz and ALMA 470GHz data images. These were all matched with the Forbrich 2011 data and [HC2000].

'IR_matched_catalog_B7.fits'. This is the table with values from the dendrograms and fluxes from each of the original images. No longer really useful except for matching the indexes of each image ('_idx_BX') to the master ID value ('D_ID'). 

'simbad_catalog.fits'. This is the same as the above but matched with a bunch of optical/IR catalogs. 

The other tables here are probably less useful
