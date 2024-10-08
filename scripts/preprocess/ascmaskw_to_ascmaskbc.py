#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For PF preprocessing
Simple watershed mask to specific bc on edge

Created on Wed Sep 04 14:34:20 2024
@author: S.P.
"""

import argparse
import matplotlib.pyplot as plt

from parflow import Run
from parflow.tools.fs import get_absolute_path
from parflow.tools.io import load_patch_matrix_from_asc_file
from parflow.tools.builders import SolidFileBuilder

import numpy as np

### PATH TO FILES #############################################
# input watershed mask ('Unique label per basin', r.watershed output)
path_in = '/home/patras/Valmalenco/Data/DataElab/'
fw_in = "maskw.c500.v2.asc" #0,1
fw_in = path_in + fw_in

path_dem = '/home/patras/Valmalenco/Data/DataElab/'
f_dem = "hydroDEM.c500.v2.asc"
f_dem = path_dem + f_dem

# output mask
path_out = '/home/patras/Valmalenco/Data/DataElab/'
filename_out = "maskb.c500.v4.asc" #top=1, bottom=2, side=3, nodata=0
f_out = path_out + filename_out

### READ ASC MASK #####################################

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
data = np.loadtxt(fw_in, skiprows=header_rows, dtype='float64')

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

# DEM

header_rows = 6
header_info = {}
row_ite = 1
with open(f_dem, 'rt') as file_h:
     for line in file_h:
        if row_ite <= header_rows:
             line = line.split(" ", 1)
             header_info[line[0]] = float(line[1])
        else:
             break
        row_ite = row_ite+1
# read data array :: mask
dem = np.loadtxt(f_dem, skiprows=header_rows, dtype='float64')

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

# POSITION BASIN CLOSURE
#print(f'data min {dem.min()}')
minpos = np.where(data[-1,:]==1)
j_minpos = minpos[0][0]
print(j_minpos)

NY3 = 5
NX3 = 10

maxhead1 = min(dem[-NY3,j_minpos-NX3:j_minpos])
maxhead2 = max(dem[-NY3,j_minpos:j_minpos+NX3+1]) #max head that could be applied
print(f'max head that can be applied for PLT as out bc {maxhead1,maxhead2}')

fig = plt.figure(1)
ax = fig.add_subplot()
#datashape = dem.shape
# print(datashape)
#x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))

ax.plot(dem[-NY3,j_minpos-NX3:j_minpos+NX3+1])
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')
ax.set_title('section of closure bc')

plt.show()

# WRITE OUTPUT MASK #####

maskb = 2*np.ones(data.shape) # (bottom=2)
maskb[-NY3:,j_minpos-NX3:j_minpos+NX3+1] = 3 # (side=3) 
ixtop = np.where(data==1)
maskb[ixtop]=1 #*maskw[ixtop]
maskb[-1,j_minpos-NX3:j_minpos+NX3+1] = 3

fout = open(f_out,"w")

fout.write("ncols"+" "+str(ncols)+"\n")
fout.write("nrows"+" "+str(nrows)+"\n")
fout.write("cellsize"+" "+str(int(cellsize))+"\n")
fout.write("xllcorner"+" "+str(0)+"\n")
fout.write("yllcorner"+" "+str(0)+"\n")
fout.write("NODATA_value"+" "+str(0)+"\n")

#for i in range(nrows-NY3):
#    for j in range(ncols):
#        fout.write(str(int(data[i,j]))+"\n")
#for i in range(nrows-NY3,nrows):
#    for j in range(j_minpos-NX3):
#        fout.write(str(int(data[i,j]))+"\n")
#    for j in range(j_minpos-NX3,j_minpos+NX3+1):
#        if int( data[i,j])==1: #if data[i,j]
#            fout.write(str(int(data[i,j]))+"\n")
#        else:
#            fout.write(str(3)+"\n")
#    for j in range(j_minpos+NX3+1,ncols):
#        fout.write(str(int(data[i,j]))+"\n")

for i in range(nrows):
    for j in range(ncols-1):
        fout.write(str(int(maskb[i,j]))+" ")
    fout.write(str(int(maskb[i,ncols-1]))+"\n")

fout.close()
