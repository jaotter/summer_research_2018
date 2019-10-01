from dendrogram_catalog import gauss_sizes_imgs

#B7 data
dirB7 = '/lustre/aoc/students/jotter/directory/OrionB7/'
nm = 'Orion_SourceI_B7_continuum_'
end = '.image.tt0.pbcor.fits'
images = [nm+'r-2.clean0.1mJy.1000klplus.deepmask'+end, nm+'r-2.clean0.5mJy.500klplus'+end, nm+'r0.5.clean0.5mJy.150mplus.deepmask'+end, nm+'r0.5.clean0.5mJy.50mplus.deepmask'+end, nm+'r2.clean1mJy.500klplus.deepmask'+end, nm+'r2.clean1mJy.50mplus.deepmask'+end, nm+'r2.clean2mJy.allbaselines'+end]
#gauss_sizes_imgs(9, 'B7', images, dirB7)

#B6 data
dirB6 = '/lustre/aoc/students/jotter/directory/OrionB6/'
nm = 'Orion_SourceI_B6_continuum_'
end = '.image.tt0.pbcor.fits'
images = [nm+'r-2.clean0.1mJy.1000klplus.deepmask'+end, nm+'r-2.clean0.5mJy.500klplus'+end, nm+'r0.5.clean0.5mJy.150mplus.deepmask'+end, nm+'r0.5.clean0.5mJy.50mplus.deepmask'+end, nm+'r2.clean1mJy.500klplus.deepmask'+end, nm+'r2.clean1mJy.50mplus.deepmask'+end, nm+'r2.clean2mJy.allbaselines'+end]
#gauss_sizes_imgs(9, 'B6', images, dirB6)

#B3 data
directory = '/lustre/aoc/students/jotter/directory/OrionB3/'
nm = 'Orion_SourceI_B3_continuum_'
end = '.image.tt0.pbcor.fits'
images = [nm+'r-2.clean0.1mJy.1000klplus.deepmask'+end, nm+'r-2.clean0.5mJy.500klplus'+end, nm+'r0.5.clean0.5mJy.150mplus.deepmask'+end, nm+'r0.5.clean0.5mJy.allbaselines.deepmask'+end, nm+'r2.clean1mJy.500klplus.deepmask'+end, nm+'r2.clean1mJy.50mplus.deepmask'+end, nm+'r2.clean2mJy.allbaselines'+end]
gauss_sizes_imgs(9, 'B3', images, directory)
