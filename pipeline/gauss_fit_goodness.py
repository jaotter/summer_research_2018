from astropy.io import fits
from astropy.table import Table
import glob

catalogs = sorted(glob.glob('/lustre/aoc/students/jotter/dendro_catalogs/*_dendro_catalog_leaves_fixed.fits'))
names = ['340GHz', '470GHz', 'B3', 'B6', 'B7_hr', 'B7_lr']
#y = good fit, n = bad fit, m = maybe, could fit with another function
GHz340_fit = ['y','n','y','y','m','y','y','n','m','y','m','y','y','y','m','m','y','y','m','y','y','y']
GHz470_fit= ['y','y','y','y','n','y','y','y','y']
B3_fit=['y','y','y','y','y','y','y','y','y','m','y','y','y','y','y','y','y','y','y','n','y','y','y','y','m','y','y','y','y','y','y','y','m','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','y','m','y','y','y','y','y','y','y']
B6_fit = ['m','m','y','m','n','y','y','y','n','y','y','n','y','n','m','y','y','m','m','y','y','m','y','y','y','y','y','y','y','y']
B7_hr_fit = ['y','b','y','n','y','n','m','n','y','y','m','m','y','y','m','y','y','n','n']
B7_lr_fit = ['y','n','y','n','m','n','y','y','y','m','y']
fit_class = [GHz340_fit, GHz470_fit, B3_fit, B6_fit, B7_hr_fit, B7_lr_fit]

for i in range(len(catalogs)):
	data = Table(fits.getdata(catalogs[i]))
	data['fit_goodness_'+names[i]] = fit_class[i]
	data.write('/lustre/aoc/students/jotter/dendro_catalogs/'+names[i]+'_dendro_catalog_leaves_fixed.fits', format='fits', overwrite=True)
