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
from mpl_toolkits.mplot3d import Axes3D
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

###
np.set_printoptions(precision=6)
###


#### READ PFB output results

parser = argparse.ArgumentParser()
parser.add_argument('runname', type = str)
parser.add_argument('variable',type = str, default="press")
parser.add_argument('dumptime',type = str, default="00010")
args = parser.parse_args()

#path_out = "/home/patras/Valmalenco/Data/DataPF/"
#filename_out = "slopeY.c500.v2.pfb"
# path_out = "/home/patras/Lombardy/Tmp/PLT_35_Leo/" #pyyamlshTmp/"
# path_out = "/home/patras/Valmalenco/Tmp/"+args.runname+"_new/"
path_out = './'
filename_out = args.runname + ".out." + args.variable + "." + args.dumptime + ".pfb"
# filename_out = "PLT_Box.out.slope_x.pfb"
# filename_out = args.pfbfilename 
f_out = path_out + filename_out
data_to_plot = read_pfb(get_absolute_path(f_out))

datashape = data_to_plot.shape
#print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)
#print(type(data_to_plot))  # numpy array

# Apply mask where data <-3.1e+38 (ie. nodata_value)
data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)


#### READ DEM


# input .asc
path_dem = '/home/patras/PF-Valmalenco/Data/DataElab/'
filename_dem = 'hydroDEM.c250.v2.asc'
f_dem = path_dem + filename_dem

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
# read data array :: dem
dem = np.loadtxt(f_dem, skiprows=header_rows, dtype='float64')
# print(dem)

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])
cellsize = header_info['cellsize']
xwest = header_info['xllcorner'] # x lower
xright = xwest+ncols*cellsize
ysouth = header_info['yllcorner'] # y lower
ynorth = ysouth+nrows*cellsize

#### PLOT

fig = plt.figure()
ax = fig.add_subplot(111, projection=Axes3D.name)

x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
y = y.max() - y

# LW_surface_press
#DZ = 2.0
#dZScale = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.05, 0])

# PLT_box
#DZ = 2.0
#dZScale = np.array([1.0,0.0]) #2.0
#DepthCenteredCellZ = np.array([-1.0])

# EI_box
DZ = 1.0
#dZScale = np.array([1.0,1.0,0.0])
#dZScale = np.array([1.7, 0.3, 0.0])
dZScale = np.array([1.0,0.4,0.3,0.15,0.1,0.05,0])

# Lombardy
#DZ = 25.0
#dZScale = np.array([5.0, 2.8, 0.8, 0.32, 0.04, 0.028, 0.006, 0.004, 0.002,0.0])

DZScaled = DZ*dZScale
DepthFaceCellZ = np.array([sum(DZScaled[i:]) for i in range(len(DZScaled))])
DepthCenteredCellZ = (-1)*(DepthFaceCellZ[:-1]+DepthFaceCellZ[1:])/2 #0 @ surface

p = -1 #surface layer
z_agl = DepthCenteredCellZ[p]
data = data_masked[p]

ax.plot_surface(x,y,,np.full_like(data, z_agl),cmap=plt.cm.Blues, rstride=1, cstride=1, antialiased=True, shade=False)

#ax.set_ylim(0,102)
ax.view_init(40, -60, 0) # elev, azimuth, roll
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z [m asl.]')
ax.set_title('Pressure head respect to surface layer [m]')

#m = plt.cm.ScalarMappable(cmap=plt.cm.Blues)  #, norm=surf.norm)
vmin = data.min()
vmax = data.max()
#ax.set_zlim(-1+DepthCenteredCellZ[0], 1+DepthCenteredCellZ[-1])
#m.set_clim(vmin,vmax)
#m.set_clim(53.9982,53.9983)
#m.set_clim(0,50)
print(vmin,vmax)
#print(f'non all nul in cell : {np.where(data>-3.402823466385288e+38)}')
#plt.colorbar(m, ax=plt.gca())

plt.show()
plt.savefig('/mnt/c/Users/Sophie/Documents/4-Figures/'+args.runname+'.'+args.variable+'.'+args.dumptime+'.3dproj.png')
# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
