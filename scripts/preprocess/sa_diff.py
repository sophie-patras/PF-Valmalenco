#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
To compare 2 files of same dimensions
sa file with 1 line header (nrows*ncols lines), to plot as surface in 3D plan

Created on Tue Sep 03 13:34:20 2024
@author: S.P.
"""

import numpy as np
import argparse
import matplotlib.pyplot as plt

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import write_pfb, read_pfb

# arg parse

parser = argparse.ArgumentParser()
parser.add_argument('safilename1', type = str)
parser.add_argument('safilename2', type = str)
#parser.add_argument('cellsize', type = float, default=500.0)
#parser.add_argument('-parflow-directory', '--parflow-directory')
#parser.add_argument('-working-directory', '--working-directory', default="~/Valmalenco/Codes/RunProcess/Tmp")
#parser.add_argument('--show-line-error', action='store_true')
#parser.add_argument('--validation-verbose', action='store_true')
# parser.add_argument('--write-yaml',action='store_true')
#parser.add_argument('-p', '--p', type=int, default=2)
#parser.add_argument('-q', '--q', type=int, default=2)
#parser.add_argument('-r', '--r', type=int, default=1)
args = parser.parse_args()

# input .sa
# path_in1 = '/home/patras/Valmalenco/Data/DataElab/'
# path_in = '/home/patras/PF-Test/'
path_in2 = '/home/patras/Lombardy/Data/SA/'
path_in1 = path_in2

filename_in1 = args.safilename1
filename_in2 = args.safilename2

f_in1 = path_in1 + filename_in1
f_in2 = path_in2 + filename_in2

# File1
header_rows = 1
#header_info = {}
#row_ite = 1
#with open(f_in, 'rt') as file_h:
#     for line in file_h:
#        if row_ite <= header_rows:
#             line = line.split(" ", 1)
#             header_info[line[0]] = float(line[1])
#        else:
#             break
#        row_ite = row_ite+1

f = open(f_in1,'r')
data1 = np.loadtxt(f_in1, skiprows=header_rows, dtype='float64')
# print(dem)
f.close()

# File2
header_rows = 1
#header_info = {}
#row_ite = 1
#with open(f_in, 'rt') as file_h:
#     for line in file_h:
#        if row_ite <= header_rows:
#             line = line.split(" ", 1)
#             header_info[line[0]] = float(line[1])
#        else:
#             break
#        row_ite = row_ite+1

f = open(f_in2,'r')
data2 = np.loadtxt(f_in2, skiprows=header_rows, dtype='float64')
f.close()

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = 265 #48 91 # int(header_info['ncols'])
nrows = 165 #49 70 # int(header_info['nrows'])
cellsize = 2000 #500 2000 #90 # header_info['cellsize']
xwest = 0 # header_info['xllcorner'] # x lower left
xright = xwest+ncols*cellsize
ysouth = 0 # header_info['yllcorner'] # y lower left
ynorth = ysouth + nrows*cellsize # header_info['yllcorner']+nrows*cellsize

map_extent = (xwest, xright, ysouth, ynorth)
print(f'map extent = {map_extent}')
#nodata_value = header_info['NODATA_value']

# data:: diff 1.sa - 2.sa
data = data1 - data2

# either correct here .sa orientation or in the y definition for plot
datanp = np.empty((nrows,ncols))
for i in range(nrows):
    for j in range(ncols):
        datanp[i,j] = data[(ncols*i)+j]
print(datanp)
data = datanp

## Plot Dem : error in colormap edges, to be corrected (how ?)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

datashape = data.shape
print(datashape)

x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))
# print(x)
x = xwest + x*cellsize
# print(x)

# print(y)
y = ysouth + y*cellsize
# print(y)

surf = ax.plot_surface(x,y,data,cmap=plt.cm.viridis, rstride=1, cstride=1, antialiased=False, shade=False)
#facecolors True
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
#ax.set_zlabel('id')
#ax.set_title('Mask')                                                   

m = plt.cm.ScalarMappable(cmap=plt.cm.viridis) #, norm=surfy.norm)
vmin = data.min()
vmax = data.max()
print(f'min,max diff : {vmin,vmax}')
ax.set_zlim(vmin,vmax)
m.set_clim(vmin,vmax)
plt.colorbar(m, ax=plt.gca())
#fig.colorbar(surf, shrink=0.5, aspect=7) #ax=plt.gca())

plt.show()
