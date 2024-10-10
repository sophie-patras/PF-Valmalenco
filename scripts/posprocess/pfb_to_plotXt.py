#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
Temporal evolution of chosen variable, at control points (defined in controlpoint.csv in EPSG:32632)
    - read asc dem (in EPSG:32632)
    - location of the control point in (i,j) coordinates
    - extra in for loop on dumptimes, variable on all layers
    - plot : 1) variable(t) for all layers for all control points (3 subplots)
             2) variable(t) of surface layer for all control points (1 plot)

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""

# import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

###
np.set_printoptions(precision=6)
FIG_DIR = '/mnt/c/Users/Sophie/Documents/4-Figures/'
plt.style.use('~/PF-Valmalenco/scripts/config.mplstyle')
###


#### READ PFB output results

parser = argparse.ArgumentParser()
parser.add_argument('runname', type = str)
parser.add_argument('variable',type = str, default="press")
parser.add_argument('dtimes',type=int,default=7)
args = parser.parse_args()

dumptimes = args.dtimes

# Input files
# .pfb files path
path_out = './'
# input .asc dem
path_dem = '/home/patras/PF-Valmalenco/data/prepareddata/'
filename_dem = 'hydroDEM.c250.v2.asc'
f_dem = path_dem + filename_dem
# controlpoints
f_cp = '/home/patras/PF-Valmalenco/data/controlpoints.txt'


# READ FILES (dem, cp, pfb)

# dem
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
# dem = np.loadtxt(f_dem, skiprows=header_rows, dtype='float64')
# print(dem)

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])
cellsize = header_info['cellsize']
xwest = header_info['xllcorner'] # x lower
xright = xwest+ncols*cellsize
ysouth = header_info['yllcorner'] # y lower
ynorth = ysouth+nrows*cellsize

# read cp

data_cp = pd.read_csv(f_cp)
x_cp = np.array(data_cp['xcoord'])
y_cp = np.array(data_cp['ycoord'])
# pfb read from lower left
i_cp = [int((y_cp[k] - ysouth)/cellsize) for k in range(3)]
j_cp = [int((x_cp[k] - xwest)/cellsize) for k in range(3)]

### DUMPTIME Loop

layers = 6

X10 = np.empty((dumptimes+1,layers))
X20 = np.empty((dumptimes+1,layers))
X30 = np.empty((dumptimes+1,layers))

for t in range(dumptimes+1):

    if t<=9:
        tstr = '0000'+str(t)
    elif t<=99:
        tstr = '000'+str(t)
    elif t<=999:
        tstr = '00'+str(t)

    filename_out = args.runname + ".out." + args.variable + "." + tstr + ".pfb"
    f_out = path_out + filename_out
    data = read_pfb(get_absolute_path(f_out))

    #datashape = data_to_plot.shape # NZ,NY,NX
    X10[t,:] = data[:,j_cp[0],i_cp[0]]
    X20[t,:] = data[:,j_cp[1],i_cp[1]]
    X30[t,:] = data[:,j_cp[2],i_cp[2]]

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

z_agl = DepthCenteredCellZ

if args.variable=='press':
    for t in range(dumptimes+1):
        X10[t,:] = X10[t,:] + z_agl[:]

#### PLOT

loc_cp = data_cp['localita']
colors_cp = ['r','g','b']

tarray = np.linspace(0,dumptimes,dumptimes+1)

fig = plt.figure(1)

ax = fig.add_subplot(111)

ax.plot(tarray,X10[:,-1],label=loc_cp[0])
ax.plot(tarray,X20[:,-1],label=loc_cp[1])
ax.plot(tarray,X30[:,-1],label=loc_cp[2])
ax.grid(True)
ax.legend()
ax.set_xlabel('t [dumptimes]')
ax.set_ylabel(args.variable)

#plt.show()
plt.savefig(FIG_DIR + args.runname + '.' + args.variable + '.L1evo.png')

figs, axs = plt.subplots(1,3,figsize=(15,5))

for l in range(6):
    axs[0].plot(tarray,X10[:,l],label='layer '+str(layers-l))
    axs[1].plot(tarray,X10[:,l],label='layer '+str(layers-l))
    axs[2].plot(tarray,X10[:,l],label='layer '+str(layers-l))

axs[1].set_xlabel('t [dumptimes]')
axs[0].set_ylabel(args.variable)
axs[0].set_title(loc_cp[0])
axs[1].set_title(loc_cp[1])
axs[2].set_title(loc_cp[2])
axs[2].legend()

plt.savefig(FIG_DIR + args.runname + '.' + args.variable + '.LAevo.png')
