#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For PF preprocessing
Slope Calculation 
    - percent
    - upwind scheme for OverlandKinematic BC

Created on Fri Jun 21 10:02:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import write_pfb, read_pfb

# input DEM
path_in = '/home/patras/Valmalenco/Data/DataElab/'
filename_in = 'hydroDEM.c250.v2.asc'
#path_in = '/home/patras/Lombardy/Data/TXT/'
#filename_in = 'DEM_35.txt'
f_in = path_in + filename_in
fx_out = path_in + 'slopeX.c250.v1.asc'
fy_out = path_in + 'slopeY.c250.v1.asc'

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

xright = xwest+ncols*cellsize
ynorth = ysouth+nrows*cellsize
map_extent = (xwest, xright, ysouth, ynorth)

### SLOPE
slopex = np.empty((nrows,ncols),dtype=float)
slopey = np.empty((nrows,ncols),dtype=float)

# upwind scheme, with starting point bottom left
# Note:
# python values are ordered from top left
# (parflow start calculating from bottom left)

# Internal cells
for i in range(nrows):
    slopex[i,:-1] = (dem[i,1:]-dem[i,:-1])/cellsize
for j in range(ncols):
    slopey[1:,j] = (dem[:-1,j]-dem[1:,j])/cellsize

# right and top border cells - identical to neighbour
#slopex[:,:-1] = (dem[:,1:]-dem[:,:-1])/cellsize
slopex[:,-1] = slopex[:,-2]
#slopey[1:,:] = (dem[:-1,:]-dem[1:,i])/cellsize
slopey[0,:] = slopey[1,:]

### WRITE

fout = open(fx_out,"w")

fout.write("ncols"+" "+str(ncols)+"\n")
fout.write("nrows"+" "+str(nrows)+"\n")
fout.write("cellsize"+" "+str(int(cellsize))+"\n")
fout.write("xllcorner"+" "+str(0)+"\n")
fout.write("yllcorner"+" "+str(0)+"\n")
fout.write("NODATA_value"+" "+str(int(nodata_value))+"\n")

for i in range(nrows):
    for j in range(ncols-1):
        fout.write(str(slopex[i,j])+" ")
    fout.write(str(slopex[i,ncols-1])+"\n")
fout.close()

#y
fout = open(fy_out,"w")

fout.write("ncols"+" "+str(ncols)+"\n")
fout.write("nrows"+" "+str(nrows)+"\n")
fout.write("cellsize"+" "+str(int(cellsize))+"\n")
fout.write("xllcorner"+" "+str(0)+"\n")
fout.write("yllcorner"+" "+str(0)+"\n")
fout.write("NODATA_value"+" "+str(int(nodata_value))+"\n")

for i in range(nrows):
    for j in range(ncols-1):
        fout.write(str(slopey[i,j])+" ")
    fout.write(str(slopey[i,ncols-1])+"\n")
fout.close()


