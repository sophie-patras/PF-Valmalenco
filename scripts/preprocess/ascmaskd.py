#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For PF preprocessing
Raster of ones, to mask the domain (maskd)
based on params given manually as input

Created on Sun Sep 8 16:34:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt

# .asc out
path_out = '/home/patras/Valmalenco/Data/DataElab/'
filename_out = 'maskd.c500.v2.asc'
f_out = path_out + filename_out

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = 48
nrows = 49
cellsize = 500.0
xwest = 0
xright = xwest+ncols*cellsize
ysouth = 0
ynorth = ysouth+nrows*cellsize
map_extent = (xwest, xright, ysouth, ynorth)

nodata_value = -9999

#data = np.empty((nrows,ncols))

# WRITE

fout = open(f_out,"w")

fout.write("ncols"+" "+str(ncols)+"\n")
fout.write("nrows"+" "+str(nrows)+"\n")
fout.write("cellsize"+" "+str(int(cellsize))+"\n")
fout.write("xllcorner"+" "+str(xwest)+"\n")
fout.write("yllcorner"+" "+str(ysouth)+"\n")
fout.write("NODATA_value"+" "+str(nodata_value)+"\n")

for i in range(nrows):
    for j in range(ncols-1):
        fout.write("1 ")#str(int(data[i,j]))+" ")
    fout.write("1\n")#str(int(data[i,ncols-1]))+"\n")
fout.close()
