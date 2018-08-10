from snr_rejection import reject_sources
from astropy.table import Table
import glob

cat_dir = '/lustre/aoc/students/jotter/dendro_catalogs/'
#catalogs = sorted(glob.glob('/lustre/aoc/students/jotter/dendro_catalogs/*_leaves.fits'))
catalogs = [cat_dir+'340GHz_dendro_catalog_leaves.fits', cat_dir+'470GHz_dendro_catalog_leaves.fits', cat_dir+'B3_dendro_catalog_leaves.fits', cat_dir+'B6_dendro_catalog_leaves.fits', cat_dir+'B7_hr_dendro_catalog_ref_B3.fits', cat_dir+'B7_lr_dendro_catalog_ref_B3.fits']


data_dir = '/lustre/aoc/students/jotter/directory/'
datname = [data_dir+'ALMA_340GHz_Nov2017_edit.fits', data_dir+'ALMA_470GHz_Aug2015.fits', data_dir+'Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits', data_dir+'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines.image.tt0.pbcor.fits', data_dir+'Orion_SourceI_B7_continuum_r-2.mask5mJy.clean4mJy.image.tt0.pbcor.fits', data_dir+'member.uid___A001_X88e_X1dd.Orion_BNKL_source_I_sci.spw25_27_29_31.cont.I.pbcor.fits']

snrs = [5,10,5,4.5,2.5,5]
#max_size_IDs = [33,None,23,53,None,None]
names = ['340GHz', '470GHz', 'B3', 'B6', 'B7_hr', 'B7_lr']


for ind in range(len(catalogs)):
	print(names[ind])
	table = reject_sources(names[ind], catalogs[ind], datname[ind], snrs[ind])
	Table(table).write('/lustre/aoc/students/jotter/dendro_catalogs/'+names[ind]+'_dendro_catalog_leaves_fixed.fits', overwrite=True)
