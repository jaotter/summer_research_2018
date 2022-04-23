import os
from astroquery.vizier import Vizier
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import coordinates
from astropy import units as u
import regions

basepath = '/Users/adam/work/students/JustinOtter/summer_research_2018/'

def getcat(catname):
    filename = os.path.join(basepath, 'tables', catname.replace("/","_")+".ecsv")
    if os.path.exists(filename):
        return Table.read(filename)
    else:
        #result = Vizier(row_limit=100000).get_catalogs(catname)[0]
        result = Vizier(row_limit=100000).query_region(SkyCoord('5:35:14.5',
                                                                '-5:22:30.6',
                                                                unit=(u.h,
                                                                      u.deg),
                                                                frame='icrs'),
                                                       radius=5*u.arcmin,
                                                       catalog=catname)[0]
        result.write(filename)
        return result

coupviz = getcat('J/ApJS/160/319/coup')
forbviz = getcat('J/ApJ/822/93')
eisner2016_gems = getcat('J/ApJ/826/16/table1')
eisner2016_mm = getcat('J/ApJ/826/16/table2')
MLLA = getcat('J/ApJ/573/366')
meingastsci = getcat('J/A+A/587/A153/science')


forbcoords = SkyCoord([x for x in forbviz['RAJ2000']],
                      [x for x in forbviz['DEJ2000']], frame='fk5', unit=(u.h,u.deg))
forbregs = regions.Regions([regions.PointSkyRegion(coord, meta={'text': id, 'shape': 'square', 'color': 'orange'})
                            for coord,id in zip(forbcoords, forbviz['Seq'])])
forbregs.write(f'{basepath}/tables/Forbrich2016.reg', overwrite=True)

coupcoords = SkyCoord([x for x in coupviz['RAJ2000']],
                      [x for x in coupviz['DEJ2000']], frame='fk5', unit=(u.deg, u.deg))

coupregs = regions.Regions([regions.PointSkyRegion(coord, meta={'text': id, 'shape': 'diamond', 'color': 'cyan'})
                            for coord,id in zip(coupcoords, coupviz['COUP'])])
coupregs.write(f'{basepath}/tables/COUP.reg', overwrite=True)


eisner2018_mm = Table.read(f'{basepath}/tables/eisner_tbl.txt', format='ascii')

crd_eis2016mm = SkyCoord(eisner2016_mm['RAJ2000'], eisner2016_mm['DEJ2000'], frame='fk5', unit=(u.h, u.deg))
crd_eis2016gems = SkyCoord(eisner2016_gems['RAJ2000'], eisner2016_gems['DEJ2000'], frame='fk5', unit=(u.h, u.deg))
crd_eis2018mm = SkyCoord(eisner2018_mm['alpha'], eisner2018_mm['delta'], frame='fk5', unit=(u.h, u.deg))

tab = B3_table = Table.read(f'{basepath}/final_tables/datafile4.txt', format='ascii.cds')
mm_coords = SkyCoord((B3_table['RAs'].quantity + (5*u.hour).to(u.s) + (35*u.min).to(u.s))*(15*u.arcsec/u.s),
                  -5*u.deg - B3_table['DEm'].quantity - B3_table['DEs'].quantity, frame='icrs', )

idx, sep, _ = crd_eis2018mm.match_to_catalog_sky(mm_coords)
new_in_e2018 = sep > 1*u.arcsec

full_mm_sourcelist = coordinates.concatenate([mm_coords.fk5, crd_eis2018mm[new_in_e2018]])
full_mm_sourcelist = SkyCoord(full_mm_sourcelist.ra, full_mm_sourcelist.dec, frame='fk5')

idx, sep, _ = crd_eis2016mm.match_to_catalog_sky(full_mm_sourcelist)
new_in_e2016 = sep > 1*u.arcsec

full_mm_sourcelist = coordinates.concatenate([full_mm_sourcelist.fk5,
                                              crd_eis2016mm[new_in_e2016].fk5])
full_mm_sourcelist = SkyCoord(full_mm_sourcelist.ra, full_mm_sourcelist.dec, frame='fk5')



idx, sep, _ = forbcoords.match_to_catalog_sky(full_mm_sourcelist)
new_in_forb = sep > 1*u.arcsec
full_radmm_sourcelist = coordinates.concatenate([full_mm_sourcelist.fk5,
                                                 forbcoords[new_in_forb].fk5])
full_radmm_sourcelist = SkyCoord(full_radmm_sourcelist.ra, full_radmm_sourcelist.dec, frame='fk5')



idx, sep, _ = coupcoords.match_to_catalog_sky(full_mm_sourcelist)
new_in_coup = sep > 1*u.arcsec
full_xradmm_sourcelist = coordinates.concatenate([full_radmm_sourcelist.fk5,
                                                  coupcoords[new_in_coup].fk5])
full_xradmm_sourcelist = SkyCoord(full_xradmm_sourcelist.ra, full_xradmm_sourcelist.dec, frame='fk5')



MLLAcoords = SkyCoord(MLLA['RAJ2000'], MLLA['DEJ2000'], frame='fk5', unit=(u.h, u.deg))
meingastcoords = SkyCoord(meingastsci['RAJ2000'], meingastsci['DEJ2000'], frame='fk5', unit=(u.h, u.deg))

idx, sep, _ = crd_eis2016gems.match_to_catalog_sky(MLLAcoords)
new_in_e2016gems = sep > 1*u.arcsec
IR_coords = coordinates.concatenate([MLLAcoords.fk5,
                                     crd_eis2016gems[new_in_e2016gems].fk5])
IR_coords = SkyCoord(IR_coords.ra, IR_coords.dec, frame='fk5')


idx, sep, _ = meingastcoords.match_to_catalog_sky(MLLAcoords)
new_in_meingast = sep > 1*u.arcsec
IR_coords = coordinates.concatenate([IR_coords.fk5,
                                     meingastcoords[new_in_meingast].fk5])
IR_coords = SkyCoord(IR_coords.ra, IR_coords.dec, frame='fk5')



idx, sep, _ = IR_coords.match_to_catalog_sky(full_mm_sourcelist)
new_in_IR = sep > 1*u.arcsec
full_irxradmm_sourcelist = coordinates.concatenate([full_xradmm_sourcelist.fk5,
                                                  IR_coords[new_in_IR].fk5])
full_irxradmm_sourcelist = SkyCoord(full_irxradmm_sourcelist.ra,
                                    full_irxradmm_sourcelist.dec, frame='fk5')

