#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical and numerical output
Pfb file (3d (p-layers)) velocity (x,y,z directions) to discharge on a cross-section, evolution. 

Application, inflow-outflow total balance check:
 - 4 plots of Q(t) of contribution of each layer + sum of all, in {x-lower,x-upper,y-lower,y-upper}
 - To evaluate the total storage vs recharge print Q_outflow/Q_recharge where
    Q_outflow = sum(Q_outtot_i) for i in {x-lower,x-upper,y-lower,y-upper}
    Q_recharge = -0.0001m * (NX*DX)*(NY*DY)

Require input : runname. 
Cellsize sets in script to 250.0 by default.
DumpTimes set in script to 31 by default.

NB: For now working only for Box domains !

Created on Mon Sep 23 12:11:23 2024
@author: S.P.
"""

# import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

###
np.set_printoptions(precision=6)


### INPUT
parser = argparse.ArgumentParser()
parser.add_argument('runname', type = str)
parser.add_argument('dumpt', type = int)
#parser.add_argument('variable',type = str, default="press")
#parser.add_argument('dumptime',type = str, default="00010")
args = parser.parse_args()

cellsize = 1000.0 #DX,DY

# LW
DZ = 2.0
dZScale = np.array([1.0,1.0,1.0,1.0,1.0,0.05])

# VM Box
#DZ = 1.0
#dZScale = np.array([1.0,1.0])
#dZScale = np.array([1.7, 0.3])
#dZScale = np.array([1.0,0.4,0.3,0.15,0.1,0.05])

# Lombardy
#DZ = 25.0
#dZScale = np.array([5.0, 2.8, 0.8, 0.32, 0.04, 0.028, 0.006, 0.004, 0.002,0.0])

# Face cells areas by layer
AScale = dZScale * cellsize

dumptimes = args.dumpt #31

###

path_in = './'
X = ['x','y']
NZ = len(dZScale)

# NB. About velocity files - values are face centered :
# velx shape : (6, 102, 99) !
# vely shape : (6, 103, 98) !

Q_xlower = np.empty((dumptimes+1,NZ+1))
Q_xupper = np.empty((dumptimes+1,NZ+1))
Q_ylower = np.empty((dumptimes+1,NZ+1))
Q_yupper = np.empty((dumptimes+1,NZ+1))

#Q_outbc = np.concatenate((Q_xlower,Q_xupper,Q_ylower,Q_yupper),axis=2)

for t in range(dumptimes+1):

    if t<=9:
        tstr = '.0000'+str(t)
    else:
        tstr = '.000'+str(t)

    for valdir in ['x']:
        
        filename_in = args.runname + ".out.vel" + valdir + tstr + ".pfb" 
        f_in = path_in + filename_in
        data_to_plot = read_pfb(get_absolute_path(f_in))
        # datashape = data_to_plot.shape
        # print(datashape). NB : (6, 102, 99)

        # Not of use for Box domain
        # Apply mask where data <-3.1e+38 (ie. nodata_value)
        data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)

        Q_xlower[t,:-1] = [sum(data_masked[p,:,0])*AScale[p] for p in range(NZ)]
        Q_xupper[t,:-1] = [sum(data_masked[p,:,-1])*AScale[p] for p in range(NZ)]
        
        Q_xlower[t,-1] = sum(Q_xlower[t,:-1])
        Q_xupper[t,-1] = sum(Q_xupper[t,:-1])

    for valdir in ['y']:

        filename_in = args.runname + ".out.vel" + valdir + tstr + ".pfb"
        f_in = path_in + filename_in
        data_to_plot = read_pfb(get_absolute_path(f_in))
        # datashape = data_to_plot.shape
        # print(datashape). NB : (6, 103, 98)

        data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)

        Q_ylower[t,:-1] = [sum(data_masked[p,0,:])*AScale[p] for p in range(NZ)]
        Q_yupper[t,:-1] = [sum(data_masked[p,-1,:])*AScale[p] for p in range(NZ)]

        Q_ylower[t,-1] = sum(Q_ylower[t,:-1])
        Q_yupper[t,-1] = sum(Q_yupper[t,:-1])

# shape of last vely file
datashape = data_to_plot.shape
ncols = datashape[2]
nrows = datashape[1]

#Last time step.
fig, ax = plt.subplots(figsize=(10.0,7.0))
x = np.linspace(1,ncols,ncols)
for n in range(NZ):
    ax.plot(x,data_masked[n,0,:], label='layer'+str(n)) # x-lower    
ax.set_xlabel('nx [/]')
ax.set_ylabel('vy [m/h]')
ax.legend()
ax.set_title('vy velocity along y-lower')

fig, ax = plt.subplots(figsize=(10.0,7.0))
y = np.linspace(1,nrows,nrows)
for n in range(NZ):
    ax.plot(y,data_masked[n,:,0], label='layer'+str(n)) # y-lower
ax.set_xlabel('ny [/]')
ax.set_ylabel('vy [m/h]')
ax.legend()
ax.set_title('vy velocity along x-lower')


Q_outbc = np.array([Q_xlower,Q_xupper,Q_ylower,Q_yupper])
Q_outlabels = ['Q_xlower', 'Q_xupper', 'Q_ylower', 'Q_yupper']
date = np.linspace(0,dumptimes,dumptimes+1)

# PRINTS
Q_recharge = 0.0001*(250.0*98)*(250.0*102)
print('Q recharge: ', Q_recharge)
Q_outflow = [sum(Q_outbc[:,t,-1]) for t in range(dumptimes)]
print('Q outflow: ', Q_outflow)

# PLOTS
side = 0
for i in range(4):
    Q = Q_outbc[i]
    fig, ax = plt.subplots(figsize=(10.0,7.0))
    for n in range(NZ):
        ax.plot(date, Q[:,n], label='layer'+str(n),marker='+',linestyle='-')
    
    ax.set_xlabel('t [h]')
    ax.set_ylabel('Q [m^3/h]')
    ax.legend()
    ax.set_title(Q_outlabels[i])
    side +=1

plt.show()
# plt.savefig('/mnt//UsersUser/Documents/POLIMI/0_TESI/')
