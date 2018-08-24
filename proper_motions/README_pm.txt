'calc_dates.py' is a few helper functions used in multiple other files, necessary in this directory.

The files using other external catalogs all rely on the output of 'match_other_catalogs.py', which is an ugly file that catalog matches a bunch of different radio/IR/optical catalogs for use in calculating proper motions. The final catalog is called 'simbad_catalog.fits'. It also creates another fits table, 'obs_dates_errs.fits'. This table contains the ra/dec errors, dates, and date error for each dataset. Requires 'calc_dates.py'.

All files named 'B3_src*' create the original plots/measurements for different sources, now outdated.

'PM_fit.py' is a helper function which takes multiple position detections and their dates to fit a proper motion direction, and then magnitude. 

'calc_proper_motions.py' was the inital proper motion calculation using Forbrich and B3/B6 data. It calculates these proper motions, and outputs a figure with proper motion vectors for sources with a pm detection. This file is outdated and not generalized at all. 

'pm_calc.py' is currently in progress and attempts to take the desired catalogs as input and output proper motion calculations for sources with detections in all the given catalogs. It seems to work for B3/B6 and Forbrich data (which is the bulk of the proper motion measurements) but not with other catalogs - needs investigation. The goal here was to create a generalized function to calculate the proper motions using a consistent method. 


