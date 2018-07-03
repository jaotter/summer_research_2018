from snr_rejection import reject_sources
from astropy.table import Table
import glob

catalogs = sorted(glob.glob('/lustre/aoc/students/jotter/dendro_catalogs/*_all.fits'))

data_dir = '/lustre/aoc/students/jotter/directory/'
datname = [data_dir+'ALMA_340GHz_Nov2017_edit.fits', data_dir+'ALMA_470GHz_Aug2015.fits', data_dir+'Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits', data_dir+'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines.image.tt0.pbcor.fits', data_dir+'member.uid___A001_X88e_X1dd.Orion_BNKL_source_I_sci.spw25_27_29_31.cont.I.pbcor.fits']

snrs = [5,10,5,4.5,2]
max_size_IDs = [39,None,53,23,15]
names = ['340GHz', '470GHz', 'B3', 'B6', 'B7']

for ind in range(len(catalogs)):
	print(names[ind])
	table = reject_sources(catalogs[ind], datname[ind], snrs[ind], max_size_IDs[ind])
	Table(table).write('/lustre/aoc/students/jotter/dendro_catalogs/'+names[ind]+'_dendro_catalog.fits', overwrite=True)
