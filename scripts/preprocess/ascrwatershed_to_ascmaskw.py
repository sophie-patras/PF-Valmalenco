#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For PF preprocessing
Simple one/no-data mask from r.watershed output (unique label mask) to asc file in list format

Created on Tue Sep 03 13:34:20 2024
@author: S.P.
"""

import argparse

from parflow import Run
from parflow.tools.fs import get_absolute_path
from parflow.tools.io import load_patch_matrix_from_asc_file
from parflow.tools.builders import SolidFileBuilder

import numpy as np

### PATH TO FILES #############################################
# input watershed mask ('Unique label per basin', r.watershed output)
path_in = '/home/patras/Valmalenco/Data/DataElab/'
fw_in = 'hydroDEM_c500v2_EPSG32632_UniqueLabPerBasin.asc'
fw_in = path_in + fw_in

# output mask
path_out = '/home/patras/Valmalenco/Data/DataElab/'
filename_out = "maskw.c500.v2.asc"
#filename_out_dummy = "maskd.c500.v1.asc"
f_out = path_out + filename_out
#f_out_d = path_out + filename_out_dummy

### READ ASC #####################################
header_rows = 6
header_info = {}
row_ite = 1
with open(fw_in, 'rt') as file_h:
     for line in file_h:
        if row_ite <= header_rows:
             line = line.split(" ", 1)
             header_info[line[0]] = float(line[1])
        else:
             break
        row_ite = row_ite+1
# read data array :: mask
watershed = np.loadtxt(fw_in, skiprows=header_rows, dtype='float64')

# # EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])
cellsize = header_info['cellsize']
xwest = header_info['xllcorner']
xeast = header_info['xllcorner']+ncols*cellsize
ysouth = header_info['yllcorner']
ynorth = header_info['yllcorner']+nrows*cellsize
print (f'map_extent = ({xwest, xeast, ysouth, ynorth}')
nodata_value = header_info['NODATA_value']

# wmask distribute values i,j growing from west,north corner
wmask = np.zeros((nrows, ncols)) + 1.*(watershed!=nodata_value) #rather *2 ?
# print(str(wmask[10,:]))

## WRITE OUTPUT MASK #####

fout = open(f_out,"w")
#fout_d = open(f_out_d, "w")

fout.write("ncols"+" "+str(ncols)+"\n")
fout.write("nrows"+" "+str(nrows)+"\n")
fout.write("cellsize"+" "+str(int(cellsize))+"\n")
fout.write("xllcorner"+" "+str(0)+"\n")
fout.write("yllcorner"+" "+str(0)+"\n")
fout.write("NODATA_value"+" "+str(0)+"\n")

#fout_d.write("ncols"+" "+str(ncols)+"\n")
#fout_d.write("nrows"+" "+str(nrows)+"\n")
#fout_d.write("cellsize"+" "+str(int(cellsize))+"\n")
#fout_d.write("xllcorner"+" "+str(0)+"\n")
#fout_d.write("yllcorner"+" "+str(0)+"\n")
#fout_d.write("NODATA_value"+" "+str(0)+"\n")

for i in range(nrows):
    for j in range(ncols-1):
        fout.write(str(int(wmask[i,j]))+" ")
    fout.write(str(int(wmask[i,ncols-1]))+"\n")
#        fout_d.write(str(1)+"\n")
fout.close()
#fout_d.close()
