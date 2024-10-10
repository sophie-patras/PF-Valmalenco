#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
Pfb file (2d (1-layer) or 3d (p-layers)) to surface layer plot

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""

# import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

###
plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#np.set_printoptions(precision=6)
###

parser = argparse.ArgumentParser()
parser.add_argument('runname', type = str)
parser.add_argument('variable',type = str, default="press")
parser.add_argument('dumptime',type = str, default="00007")
args = parser.parse_args()

variable = args.variable

path_in = './' # execute from input file folder
#filename_out = args.runname + ".out." + args.variable + "." + args.dumptime + ".pfb"
#f_out = path_out + filename_out

if variable == 'vel':
    filename_x = args.runname + ".out." + variable + "x." + args.dumptime + ".pfb"
    filename_y = args.runname + ".out." + variable + "y." + args.dumptime + ".pfb"
    filename_z = args.runname + ".out." + variable + "z." + args.dumptime + ".pfb"
    datax = read_pfb(get_absolute_path(path_in + filename_x))
    datay = read_pfb(get_absolute_path(path_in + filename_y))
    dataz = read_pfb(get_absolute_path(path_in + filename_z))
    data = np.sqrt(np.square(datax[:,:,:-1]) + np.square(datay[:,1:,:]) + np.square(dataz[:-1,:,:]))
else:
    filename_in = args.runname + ".out." + args.variable + "." + args.dumptime + ".pfb"
    f_in = path_in + filename_in
    data = read_pfb(get_absolute_path(f_in))

datashape = data.shape
#print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)

# Apply mask where data <-3.1e+38 (ie. nodata_value)
#data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)

# LW_surface_press
#DZ = 2.0
#dZScale = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.05, 0])

# PLT_box
#DZ = 2.0
#dZScale = np.array([1.0,0.0]) #2.0
#DepthCenteredCellZ = np.array([-1.0])

# Valmalenco
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

p = -1
z_agl = DepthCenteredCellZ[p] 
data = data[p]
#data = data_masked[p]

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
    varlabel = 'h [m agl]'
    data = data + z_agl
    varrange = [-2, 1]
elif variable == 'satur':
    colorofmap = 'Blues'
    varlabel = 'S [-]'
    varrange = [0, 1]

vmin = data.min()
vmax = data.max()
vmean = np.mean(data)

fig = plt.figure(dpi=150)
ax = fig.add_subplot(111)

plt.imshow(data,cmap = colorofmap, interpolation='nearest') 
plt.colorbar(label = varlabel)
#plt.clim(varrange[0],varrange[1])

#ax.set_xlabel('nx')
#ax.set_ylabel('ny')
ax.set_ylim(0, 102)

#m = plt.cm.ScalarMappable(cmap=plt.cm.Blues)  #, norm=surf.norm)
#vmin = data.min()
#vmax = data.max()

#m.set_clim(dmin,dmax)
#m.set_clim(53.9982,53.9983)
#m.set_clim(0,50)
#print(vmin,vmax)
#print(f'non all nul in cell : {np.where(data>-3.402823466385288e+38)}')
#plt.colorbar(m, ax=plt.gca())

plt.show()
plt.savefig('/mnt/c/Users/Sophie/Documents/4-Figures/'+args.runname+'.'+args.variable+'.'+args.dumptime+'.2dxy.png')
# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
