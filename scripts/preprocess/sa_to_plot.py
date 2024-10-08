#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
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
parser.add_argument('safilename', type = str)
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
path_in = '/home/patras/PF-Valmalenco/Data/DataElab/'
# path_in = '/home/patras/PF-Test/'
# path_in = '/home/patras/Lombardy/Data/SA/'
filename_in = args.safilename
f_in = path_in + filename_in

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

f = open(f_in,'r')                                                                                                                   
dem = np.loadtxt(f_in, skiprows=header_rows, dtype='float64')
print(dem)

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = 48 #48, 265, 91 # int(header_info['ncols'])
nrows = 49 #49, 165, 70 # int(header_info['nrows'])
cellsize = 500 #500, 2000, 90 # header_info['cellsize']
xwest = 0 # header_info['xllcorner'] # x lower left
xright = xwest+ncols*cellsize
ysouth = 0 # header_info['yllcorner'] # y lower left
ynorth = ysouth + nrows*cellsize # header_info['yllcorner']+nrows*cellsize

map_extent = (xwest, xright, ysouth, ynorth)
print(f'map extent = {map_extent}')
#nodata_value = header_info['NODATA_value']

# zsouthwest = dem[-1,0]
zlower = dem.min()
zupper = dem.max()
print(f'min,max = {zlower},{zupper}')

# either correct here .sa orientation or in the y definition for plot
demnp = np.empty((nrows,ncols))
for i in range(nrows):
    for j in range(ncols):
        demnp[i,j] = dem[(ncols*i)+j]
print(demnp)
dem = demnp
print(dem[0,0])

## Plot Dem : error in colormap edges, to be corrected (how ?)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

datashape = dem.shape
print(f'shape : {datashape}')

x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))
# print(x)
x = xwest + x*cellsize
# print(x)

# print(y)
y = ysouth + y*cellsize
# print(y)

surf = ax.plot_surface(x,y,dem,cmap=plt.cm.viridis, rstride=1, cstride=1, antialiased=False, shade=False)
#facecolors True
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
#ax.set_zlabel('id')
#ax.set_title('Mask')                                                   

m = plt.cm.ScalarMappable(cmap=plt.cm.viridis) #, norm=surfy.norm)
vmin = dem.min()
vmax = dem.max()
ax.set_zlim(vmin,vmax)
m.set_clim(vmin,vmax)
plt.colorbar(m, ax=plt.gca())
#fig.colorbar(surf, shrink=0.5, aspect=7) #ax=plt.gca())

plt.show()
plt.savefig('/mnt/c/Users/Sophie/Documents/4-Figures/'+args.safilename[:-3]+'.png')
