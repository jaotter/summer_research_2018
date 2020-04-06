from region_file import compute_regions
import regions
import numpy as np


def run_astrodendro(name, min_val, min_del, n_pix, img_path):
        #reg_fname = f'/home/jotter/nrao/images/dendro_regions/{name}_minval{min_val}_mindel{min_del}_npix{n_pix}.reg'
        reg_fname = f'/home/otter/otter/nrao/dendro_regions/{name}_minval{min_val}_mindel{min_del}_npix{n_pix}.reg'
        
        dendro, cat = compute_regions(min_val, min_del, n_pix, img_path, reg_fname)
        cat.rename_column('_idx', '_idx_'+name)

def test_params_astrodendro(name, min_val_arr, min_del_arr, n_pix, img_path):
        for min_val in min_val_arr:
                for min_del in min_del_arr:
                        run_astrodendro(name, min_val, min_del, n_pix, img_path)


name = 'B3'
min_val_arr = [0.00002,0.00003,0.00004,0.00005,0.00006]
min_del_arr = [0.00005,0.00007,0.00009]
n_pix = 10
img_path = '/home/otter/otter/nrao/Orion_SourceI_B3_continuum_r0.5.clean0.05mJy.allbaselines.deepmask.image.tt0.pbcor.fits'

test_params_astrodendro(name, min_val_arr, min_del_arr, n_pix, img_path)
