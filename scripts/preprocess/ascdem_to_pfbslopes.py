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

# input DEM
path_in = '/home/patras/Valmalenco/Data/DataElab/'
filename_in = 'hydroDEM.c500.v2.asc'
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
xright = header_info['xllcorner']+ncols*cellsize
ysouth = header_info['yllcorner']
ynorth = header_info['yllcorner']+nrows*cellsize
map_extent = (xwest, xright, ysouth, ynorth)

nodata_value = header_info['NODATA_value']

### SLOPE
slopex = np.empty((nrows,ncols))
slopey = np.empty((nrows,ncols))

# upwind scheme, with starting point bottom left
# Note:
# python values are ordered from top left
# (parflow start calculating from bottom left)

# Internal cells
for i in range(nrows):
    slopex[i,:-1] = (dem[i,1:]-dem[i,:-1])/cellsize
for j in range(ncols):
    slopey[1:,j] = (dem[:-1,j]-dem[1:,j])/cellsize

# right and top border cells - identical to neighbour
slopex[:,-1] = slopex[:,-2]
slopey[0,:] = slopey[1,:]

# Plot numpy slopes for verif

#fig = plt.figure(2)
#ax = fig.add_subplot(121, projection='3d')

#datashape = slopex.shape
#x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))
#x = xwest + x*cellsize
#y = ynorth - y*cellsize
#print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)

#ax.plot_surface(x,y,np.full_like(slopex, 0), facecolors=plt.cm.RdBu(slopex), rstride=1, cstride=1, antialiased=True, shade=False)
#ax.set_xlabel('X [m]')
#ax.set_ylabel('Y [m]')
#ax.set_zlabel('i')
#ax.set_title('slope_x')

#m = plt.cm.ScalarMappable(cmap=plt.cm.RdBu) #, norm=surfy.norm)
#ax.set_zlim(-1,1)
#m.set_clim(-1,1)
#plt.colorbar(m, aspect = 10, pad=0.1, fraction=0.05, ax=plt.gca())

#ax = fig.add_subplot(122, projection='3d')

#datashape = slopey.shape
##x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))
##x = xwest + x*cellsize
##y = ynorth - y*cellsize
#surfy = ax.plot_surface(x,y,np.full_like(slopey, 0), facecolors=plt.cm.RdBu(slopey), rstride=1, cstride=1, antialiased=True, shade=False)
#ax.set_xlabel('X [m]')
#ax.set_ylabel('Y [m]')
# ax.set_zlabel('i')
#ax.set_title('slope_y')

#m = plt.cm.ScalarMappable(cmap=plt.cm.RdBu) #, norm=surfy.norm)
#ax.set_zlim(-1,1)
#print(slopex.min(), slopex.max())
#print(slopey.min(), slopey.max())
#m.set_clim(-1,1)
#plt.colorbar(m, aspect = 10, pad=0.1, fraction=0.05, ax=plt.gca())

#plt.show()

### WRITE PFB OUTPUT SLOPE FILE

# write_pfb(file, array, p=1, q=1, r=1, x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0, dz=1.0, z_first=True, dist=True, **kwargs)

version = '3'

path_out = "/home/patras/Valmalenco/Data/DataPF/"
filename_out = "slopeX.c"+ str(int(cellsize)) + ".v" + version +".pfb"
f_out = path_out + filename_out
write_pfb(get_absolute_path(f_out), slopex, p=2, q=2, r=1, x=0, y=0, z=0, dx=cellsize, dy=cellsize, dz=1.0, z_first=True, dist=False)
slopexpfb = read_pfb(get_absolute_path(f_out))

filename_out = "slopeY.c"+ str(int(cellsize)) + ".v" + version +".pfb"
f_out = path_out + filename_out
write_pfb(get_absolute_path(f_out), slopey, p=2, q=2, r=1, x=0, y=0, z=0, dx=cellsize, dy=cellsize, dz=1.0, z_first=True, dist=False)
slopeypfb = read_pfb(get_absolute_path(f_out))

# Plot pfb slopes for verif
fig = plt.figure(3)
ax = fig.add_subplot(121, projection='3d')

datashape = slopexpfb.shape
x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
x = xwest + x*cellsize
y = ynorth - y*cellsize
#print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)

ax.plot_surface(x,y,np.full_like(slopexpfb[0], 0), facecolors=plt.cm.RdBu(slopex), rstride=1, cstride=1, antialiased=True, shade=False)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
#ax.set_zlabel('i')
ax.set_title('pfb slope_x')

m = plt.cm.ScalarMappable(cmap=plt.cm.RdBu) #, norm=surfy.norm)
ax.set_zlim(-1,1)
m.set_clim(-1,1)
plt.colorbar(m, aspect = 10, pad=0.1, fraction=0.05, ax=plt.gca())

ax = fig.add_subplot(122, projection='3d')

datashape = slopeypfb.shape
#x,y = np.meshgrid(np.arange(datashape[1]), np.arange(datashape[0]))
#x = xwest + x*cellsize
#y = ynorth - y*cellsize
ax.plot_surface(x,y,np.full_like(slopeypfb[0], 0), facecolors=plt.cm.RdBu(slopey), rstride=1, cstride=1, antialiased=True, shade=False)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
# ax.set_zlabel('i')
ax.set_title('pfb slope_y')

m = plt.cm.ScalarMappable(cmap=plt.cm.RdBu) #, norm=surfy.norm)
#ax.set_zlim(vmin,vmax)
ax.set_zlim(-1,1)
print(slopex.min(), slopex.max())
print(slopey.min(), slopey.max())
m.set_clim(-1,1)
plt.colorbar(m, aspect = 10, pad=0.1, fraction=0.05, ax=plt.gca())

plt.show()
