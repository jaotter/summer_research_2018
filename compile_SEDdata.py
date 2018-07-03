#This file compiles the various SEDdata_TxtFiles into one list of lists with the HC2000 and RRS2008 IDs
from astropy.table import Table,hstack
import glob
import numpy as np
import csv

SED_files = glob.glob('/lustre/aoc/students/jotter/SEDdata_TxtFiles/*')
comp_tab = []

for f in sorted(SED_files):
    file = open(f)
    line = file.readline()
    IDs = line.split()
    HC_ID = 'none'
    RRS_ID = 'none'
    for i in IDs:
        if '[HC2000]' in i:
            HC_ID = i[i.find('[HC2000]')+8:].strip()
        if '[RRS2008]' in i:
            RRS_ID = i[i.find('[RRS2008]')+9:].strip()
    tab = np.loadtxt(f,skiprows=1,dtype='str')
    tab = np.append([HC_ID], np.append([RRS_ID], tab.reshape(-1)))
    comp_tab.append(tab)

myfile = open('/lustre/aoc/students/jotter/SED_list.csv', 'w')
with myfile:  
   writer = csv.writer(myfile)
   writer.writerows(comp_tab)
