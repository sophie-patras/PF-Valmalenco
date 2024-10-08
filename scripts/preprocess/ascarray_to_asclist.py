#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For file format conversion .asc to .asc in list format (not very useful...)

Created on Tue Sep 03 13:34:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt

# FILES PATH
path_in = '/home/patras/Valmalenco/Data/DataElab/'
filename_in = 'hydroDEM.c500.v2.asc'
f_in = path_in + filename_in

filename_out = 'hydroDEM.c500.v2.sa'
f_out = path_in + filename_out

# READ INPUT

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
# print(dem)

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])
cellsize = header_info['cellsize']
xwest = header_info['xllcorner'] # x lower left
xright = header_info['xllcorner']+ncols*cellsize
ysouth = header_info['yllcorner'] # y lower left
ynorth = header_info['yllcorner']+nrows*cellsize

map_extent = (xwest, xright, ysouth, ynorth)
print(f'map extent = {map_extent}')
nodata_value = header_info['NODATA_value']

# WRITE FILE

fout = open(f_out,"w")

fout.write("ncols"+" "+str(ncols)+"\n")
fout.write("nrows"+" "+str(nrows)+"\n")
fout.write("cellsize"+" "+str(int(cellsize))+"\n")
fout.write("xllcorner"+" "+str(0)+"\n")
fout.write("yllcorner"+" "+str(0)+"\n")
fout.write("NODATA_value"+" "+str(0)+"\n")

for i in range(nrows):
    for j in range(ncols):
        fout.write(str(int(dem[i,j]))+"\n")
fout.close()
                     
