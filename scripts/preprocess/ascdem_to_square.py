#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For PF preprocessing
DEM rectangle (used for original DEM that was reprojected and not full lined) 

Created on Fri Sep 13 10:02:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import write_pfb, read_pfb

# input DEM
path_in = '../../data/prepareddata/'
filename_in = 'hydroDEM.c250.v1.asc'
#path_in = '/home/patras/Lombardy/Data/TXT/'
#filename_in = 'DEM_35.txt'
f_in = path_in + filename_in
f_out = path_in + 'hydroDEM.c250.v2.asc'

header_rows = 6
header_info = {}
row_ite = 1
with open(f_in, 'rt') as file_h:
     for line in file_h:
        if row_ite <= header_rows:
             line = line.split(" ", 1)
             header_info[line[0]] = float(line[1])
        else:
             break
        row_ite = row_ite+1
# read data array :: dem
dem = np.loadtxt(f_in, skiprows=header_rows, dtype='float64')

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])
cellsize = header_info['cellsize']
xwest = header_info['xllcorner']
ysouth = header_info['yllcorner']
nodata_value = header_info['NODATA_value']

# Lombardy
#ncols = 265 #int(header_info['ncols'])
#nrows = 165 #int(header_info['nrows'])
#cellsize = 2000.0 #header_info['cellsize']
#xwest = 0 #header_info['xllcorner']
#ysouth = 0  #header_info['yllcorner']
#nodata_value = -9999 #header_info['NODATA_value']

xeast = xwest+ncols*cellsize
ynorth = ysouth+nrows*cellsize
map_extent = (xwest, xeast, ysouth, ynorth)
print('original :',map_extent)
#print((xeast-xwest)/cellsize)
#print((ynorth-ysouth)/cellsize)

# manually
dembis = dem[1:-1,1:-1]

xwest = xwest + cellsize
ysouth = ysouth + cellsize
xeast = xeast - cellsize
ynorth = ynorth - cellsize
map_extent = (xwest, xeast, ysouth, ynorth)
print('resize :',map_extent)

### WRITE

fout = open(f_out,"w")

fout.write("ncols"+" "+str(ncols-2)+"\n")
fout.write("nrows"+" "+str(nrows-2)+"\n")
fout.write("cellsize"+" "+str(cellsize)+"\n")
fout.write("xllcorner"+" "+str(xwest)+"\n")
fout.write("yllcorner"+" "+str(ysouth)+"\n")
fout.write("NODATA_value"+" "+str(int(nodata_value))+"\n")

for i in range(nrows-2):
    for j in range(ncols-3):
        fout.write(str(dembis[i,j])+" ")
    fout.write(str(dembis[i,ncols-3])+"\n")
fout.close()
