#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
Pfb file (2d (1-layer) or 3d (p-layers)) to stacked plot

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""

# import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

parser = argparse.ArgumentParser()
parser.add_argument('pfbfilename', type = str)
args = parser.parse_args()

#path_out = "/home/patras/Valmalenco/Data/DataPF/"
#filename_out = "slopeY.c500.v2.pfb"
# path_out = "/home/patras/Lombardy/Tmp/PLT_35_Leo/" #pyyamlshTmp/"
path_out = "/home/patras/Valmalenco/Tmp/PLT_v1_new/"
# filename_out = "PLT_Box.out.press.00000.pfb"
# filename_out = "PLT_Box.out.slope_x.pfb"
filename_out = args.pfbfilename 
f_out = path_out + filename_out
data_to_plot = read_pfb(get_absolute_path(f_out))

datashape = data_to_plot.shape
#print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)
#print(type(data_to_plot))  # numpy array

######################################

# Dem elevation
path_in = '/home/patras/Valmalenco/Data/DataElab/'
filename_in = 'hydroDEM.c500.v2.asc'
#path_in = '/home/patras/Lombardy/Data/TXT/'
#filename_in = 'DEM_35.txt'
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

##################################

demreversed = np.empty(data_to_plot.shape)
demreversed[0] = [dem[-i,:] for i in range(dem.shape[0])]

press = data_to_plot - demreversed
print(press) 
#################################

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
#y = y.max() - y

for p in range(datashape[0]):
	data = press[p]
	ax.plot_surface(x,y,np.full_like(data, p), facecolors=plt.cm.Blues(data), rstride=1, cstride=1, antialiased=True, shade=False)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('i')
ax.set_title('pressure (m)')

m = plt.cm.ScalarMappable(cmap=plt.cm.Blues)  #, norm=surf.norm)
vmin = data.min()
vmax = data.max()
ax.set_zlim(vmin,vmax)
m.set_clim(max(-200,vmin),vmax)
print(vmin,vmax)
print(f'non all nul in cell : {np.where(data>-3.402823466385288e+38)}')
plt.colorbar(m, ax=plt.gca())

plt.show()
# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
