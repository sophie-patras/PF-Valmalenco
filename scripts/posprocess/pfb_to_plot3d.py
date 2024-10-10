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

parser = argparse.ArgumentParser()
parser.add_argument('runname', type = str)
parser.add_argument('variable',type = str, default="press")
parser.add_argument('dumptime',type = str, default="00010")
args = parser.parse_args()

variable = args.variable

# input file
path_out = './'
filename_out = args.runname + ".out." + args.variable + "." + args.dumptime + ".pfb"
f_out = path_out + filename_out

data_to_plot = read_pfb(get_absolute_path(f_out))

datashape = data_to_plot.shape
#print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)
#print(type(data_to_plot))  # numpy array

# Apply mask where data <-3.1e+38 (ie. nodata_value)
data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection=Axes3D.name)

x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
#y = y.max() - y

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

#print('NB. Only 3 layers ploted')
dl = max(1,int(datashape[0]/3))

dmin = data_masked[-1].min()
dmax = data_masked[-1].max()

if variable[:-1]=='vel':
    colorofmap = 'Spectral'
    varlabel = '$v_' + variable[-1] + '$ [m/s]'
    data = data * 1/(60*60)
    varrange = [data.min(),data.max()]
elif variable == 'vel':
    colorofmap = 'Spectral'
    varlabel = 'v [m/s]'
    data = data * 1/(60*60)
elif variable == 'press':
    colorofmap = 'Blues'
    varlabel = 'h [m agl.]'
    #data = data + z_agl
    varrange = [-2, 1]
elif variable == 'satur':
    colorofmap = 'Blues'
    varlabel = 'S [-]'
    varrange = [0, 1]


for p in range(0,datashape[0]):
    Zvssurf = DepthCenteredCellZ[p]
    if (args.variable=='press'):
        data = data_masked[p] + Zvssurf
    else:
        data = data_masked[p] 
    #print(data)
    vmin = data.min()
    vmax = data.max()
    vmean = np.mean(data)
    print(f'layer {datashape[0]-p} centered depth z=',round(Zvssurf,3),'m.agl; vmin,vmean,vmax',round(vmin,5),round(vmean,5),round(vmax,5))
    if p in range(0,datashape[0],1):
        ax.plot_surface(x,y,np.full_like(data, Zvssurf), facecolors=plt.cm.Blues(data), rstride=1, cstride=1, antialiased=True, shade=False)
    #if vmin<dmin: dmin=vmin
    #if vmax>dmax: dmax=vmax

ax.view_init(30, -70, 0) # elev, azimuth, roll
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z') # [m agl.]')
#ax.set_title('stacked '+args.variable)
#ax.text2D(0.1,0.9,f"time={DumpInt*DumpGap}h",transform = ax.transAxes)

m = plt.cm.ScalarMappable(cmap=plt.cm.Blues)  #, norm=surf.norm)
#vmin = data.min()
#vmax = data.max()
#ax.set_zlim(-1+DepthCenteredCellZ[0], 1+DepthCenteredCellZ[-1])
m.set_clim(dmin,dmax)
#m.set_clim(53.9982,53.9983)
#m.set_clim(0,50)
#print(vmin,vmax)
#print(f'non all nul in cell : {np.where(data>-3.402823466385288e+38)}')
plt.colorbar(m, ax=plt.gca(), label=varlabel)

plt.show()
plt.savefig('/mnt/c/Users/Sophie/Documents/4-Figures/'+args.runname+'.'+args.variable+'.'+args.dumptime+'.png')
# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
