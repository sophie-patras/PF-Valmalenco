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

plt.style.use('../config.mplstyle')

# path to files
path_in = "../../data/prepareddata/" # corrected EPSG:32632 and framed to hydroDem.c250 with QGIS
param = "KS" # in {KS, THS, MRC}
fname_in = "_sl4_VM.asc"
f_in = path_in + param + fname_in

path_out = "/home/patras/PF-Valmalenco/data/pfdata/"
version = str(1)

# read input

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
data = np.loadtxt(f_in, skiprows=header_rows, dtype='float64')

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])
cellsize = header_info['cellsize']
xwest = header_info['xllcorner']
xright = header_info['xllcorner']+ncols*cellsize
ysouth = header_info['yllcorner']
ynorth = header_info['yllcorner']+nrows*cellsize
map_extent = (xwest, xright, ysouth, ynorth)

nodata_value = header_info['NODATA_value']

# mask where no data (glacier and permafrost)
idxmask = np.where(data == nodata_value)

if param == 'KS':
    data = data*1/100*0.01/24 # [m/h]
    data[idxmask] = data.max()
    data[idxmask] = data.min()

print(data)
# data = np.ma.masked_where(data == nodata_value, data, copy=True)

idxmask = np.where(data == nodata_value)
print(idxmask)

#data = data[::-1,::-1] #rotate by 180 degree
data = np.flip(data,0)

# write pfb
# write_pfb(file, array, p=1, q=1, r=1, x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0, dz=1.0, z_first=True, dist=True, **kwargs)

fname_out = param + ".c" + str(int(cellsize)) + ".v" + version +".pfb"
f_out = path_out + fname_out
write_pfb(get_absolute_path(f_out), data, p=4, q=4, r=1, x=0, y=0, z=0, dx=cellsize, dy=cellsize, dz=1.0, z_first=True, dist=False)

