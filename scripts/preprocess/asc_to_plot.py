#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
asc file with 6 lines header and data in matrix format, to plot as surface in 3D plan

Created on Tue Sep 03 13:34:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import write_pfb, read_pfb

# arg parse

parser = argparse.ArgumentParser()
parser.add_argument('ascfilename', type = str)
args = parser.parse_args()

# input .asc
path_in = '/home/patras/PF-Valmalenco/Data/DataElab/'
# path_in = '/home/patras/Lombardy/Data/TXT/'
filename_in = args.ascfilename
#filename_in = 'hydroDEM.c500.v2.asc'
#filename_in = 'maskw.c500.v1.asc'
f_in = path_in + filename_in

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
xwest = header_info['xllcorner'] # x lower
xright = xwest+ncols*cellsize
ysouth = header_info['yllcorner'] # y lower
ynorth = ysouth+nrows*cellsize

map_extent = (xwest, xright, ysouth, ynorth)
print(f'map extent = {map_extent}')
nodata_value = -9999 # header_info['NODATA_value']

#zsouthwest = dem[-1,0]
zlower = dem.min()
zupper = dem.max()
print(f'dem lower,upper = {zlower},{zupper}')

## Plot Dem : error in colormap edges, to be corrected (how ?)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

datashape = dem.shape
# print(datashape)
x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))
# print(x)
x = xwest + x*cellsize
# print(x)

# print(y)
y = ynorth - y*cellsize
# print(y)

surf = ax.plot_surface(x,y,dem,cmap=plt.cm.viridis, rstride=1, cstride=1, antialiased=False, shade=False)
#facecolors True
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
#ax.set_zlabel('Z [m asl.]')
#ax.set_title('DEM')                                                   
m = plt.cm.ScalarMappable(cmap=plt.cm.viridis) #, norm=surfy.norm)
vmin = dem.min()
vmax = dem.max()
ax.set_zlim(vmin,vmax)
m.set_clim(vmin,vmax)
plt.colorbar(m, ax=plt.gca())
#fig.colorbar(surf, shrink=0.5, aspect=7) #ax=plt.gca())

plt.show()
plt.savefig('/mnt/c/Users/Sophie/Documents/4-Figures/'+args.ascfilename[:-4]+'.3d.png')
