#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
Read all times <runname>.out.<variable>.<dumptime>.pfb files (p-layers) of an assigned folder path
Functions : 
- 3d
- 2dxy
- 2dxz
- 1dz(t)
...etc.

Created on Mon Oct 28 16:41:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# import imageio

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb
from parflow import Run
from parflow.tools.hydrology import calculate_surface_storage, calculate_water_table_depth, calculate_overland_flow_grid, calculate_subsurface_storage

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'

###############################################################################################
# INPUT

#path = '/home/patras/PF-Test/DumbSquare/outputs/'
#foldername = 'DS.c1_v8'
#runname = 'DS.c1'

path = '/home/patras/PF-Test/'
foldername = 'LW'
runname = 'LW_surface_press'

###############################################################################################

run = Run.from_definition(f'{path}{foldername}/{runname}.pfidb')
data = run.data_accessor

dx = data.dx
dy = data.dy
dz = data.dz
#print(dz)

nx = data.shape[2]
ny = data.shape[1]
nz = data.shape[0]
#nt = data.time

porosity = data.computed_porosity
specific_storage = data.specific_storage
mask = data.mask
slopex = data.slope_x               # shape (ny, nx)
slopey = data.slope_y               # shape (ny, nx)
#mannings = data.mannings

def read_pfidb(path,runname):
    
    run = Run.from_definition(f'{path}{foldername}/{runname}.pfidb')
    data = run.data_accessor

    dx = data.dx
    dy = data.dy
    dz = data.dz
    #print(dz)

    nx = data.shape[2]
    ny = data.shape[1]
    nz = data.shape[0]
    nt = data.time

    porosity = data.computed_porosity
    specific_storage = data.specific_storage
    mask = data.mask
    slopex = data.slope_x               # shape (ny, nx)
    slopey = data.slope_y               # shape (ny, nx)
    mannings = data.mannings

def layers_centereddepth():
    
    dzscaled = np.concatenate((dz,[0]))
    facedepth = np.array([sum(dzscaled[i:]) for i in range(nz+1)])
    centereddepth = (-1)*(facedepth[:-1]+facedepth[1:])/2 #0 @ surface
    #print(centereddepth)
    return centereddepth

def list_dumptimes(): #(path,runname):
    ndt = 7
    return ndt

def readpfblist_to_4Darray(path,runname,variable):

    nt = list_dumptimes() #path,runname)
    variable4D_rawpfb = np.empty((nt,nz,ny,nx))
    data.time=0
    if variable == 'press':
        for t in range(nt):
            variable4D_rawpfb[t] = data.pressure
            data.time +=1
    if variable == 'satur':
        for t in range(nt):
            variable4D_rawpfb[t] = data.saturation
            data.time +=1
    if variable == 'velx': #NO !
        for t in range(nt):
            variable4D_rawpfb[t] = data.velx
            data.time +=1
    if variable == 'vely':
        for t in range(nt):
            variable4D_rawpfb[t] = data.vely
            data.time +=1
    if variable == 'velz':
        for t in range(nt):
            variable4D_rawpfb[t] = data.velz
            data.time +=1

    return variable4D_rawpfb

def plot2dxz(variable4D_rawpfb,ny,varname):

    array = variable4D_rawpfb[:,:,ny,:]
    nt = array.shape[0]
    #print(array.shape)

    x = np.linspace(0,nx-1,nx)
    z = layers_centereddepth()
    z = [round(z[i],4) for i in range(nz)]

    for t in range(nt):

        fig, ax = plt.subplots()
        plt.pcolormesh(x, z, array[t], cmap='Blues')
        varlabel = 'h [m above layer]'
        varrange = [array[t].min(), array[t].max()]

        #plt.imshow(data,cmap = colorofmap, origin="lower", extent=[0,10,-1,0], aspect=2) #, interpolation='nearest'/'none')
        
        ax.set(yticks=z, yticklabels=z)
        plt.colorbar(label = varlabel)
        #plt.colorbar(orientation='vertical', label = r'pressure head vs sl1 [m]', fraction = 0.04, pad = 0.0)
        plt.clim(varrange[0],varrange[1])
        
        #plt.axis('off')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('z [m agl]')
        plt.savefig(f'{path_fig}{foldername}.{varname}.{str(t).zfill(5)}.2dxz.png', dpi = 300)
        plt.close()

def plot2dxy(variable4D_rawpfb,nz,varname):

    array = variable4D_rawpfb[:,nz,:,:]
    nt = array.shape[0]

    x = np.linspace(0,nx-1,nx)
    y = np.linspace(0,ny-1,ny)

    for t in range(nt):

        fig, ax = plt.subplots()
        plt.pcolormesh(x, y, array[t], cmap='Blues')
        varlabel = 'h [m above layer]'
        varrange = [array.min(), array.max()]

        #plt.imshow(data,cmap = colorofmap, origin="lower", extent=[0,10,-1,0], aspect=2) #, interpolation='nearest'/'none')
        
        #ax.set(yticks=z, yticklabels=z)
        plt.colorbar(label = varlabel)
        #plt.colorbar(orientation='vertical', label = r'pressure head vs sl1 [m]', fraction = 0.04, pad = 0.0)
        plt.clim(varrange[0],varrange[1])
        
        #plt.axis('off')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        plt.savefig(f'{path_fig}{foldername}.{varname}.{str(t).zfill(5)}.2dxy.png', dpi = 300)
        plt.close()

def variablescaling(array,variable):

    if variable[:-1]=='vel':
        colorofmap = 'Spectral'
        varlabel = '$v_' + variable[-1] + '$ [m/h]'
        varrange = [array.min(),array.max()]
    elif variable == 'vel':
        colorofmap = 'Spectral'
        varlabel = 'v [m/h]'
        varrange = [array.min(),array.max()]
    elif variable == 'press':
        colorofmap = 'Blues'
        varlabel = 'h [m above layer]'
        varrange = [array.min(), array.max()] #min(data.max(),2)]
        #varrange = [-0.7, 0.1]
    elif variable == 'satur':
        colorofmap = 'Blues'
        varlabel = 'S [-]'
        #varrange = [data.min(), data.max()]
        varrange = [0, 1]

    return colorofmap, varlabel, varrange

def plotmosaic(runname,projection,idx):

    if projection == '2dxz':
        x = np.linspace(0,nx-1,nx)
        y = layers_centereddepth()
        y = [round(y[i],4) for i in range(nz)]
    if projection == '2dxy':
        x = np.linspace(0,nx-1,nx)
        y = np.linspace(0,ny-1,ny)
    if projection == '2dyz':
        x = np.linspace(0,ny-1,ny)
        y = layers_centereddepth()
        y = [round(y[i],4) for i in range(nz)]

    fig, axs = plt.subplots(nt,len(pfb_outputvariables),figsize=(len(pfb_outputvariables)*5,nt*5+2))

    for v in range(len(pfb_outputvariables)):
        vname = pfb_outputvariables[v]
        array = readpfblist_to_4Darray(path,runname,vname)
        if projection == '2dxz':
            array = array[:,:,idx]
        if projection == '2dxy':
            array = array[:,idx]
        if projection == '2dyz':
            array = array[:,:,:,idx]
        pltsettings = variablescaling(array, pfb_outputvariables[v])
        for t in range(nt):
            im = axs[t,v].pcolormesh(x, y, array[t], shading='auto', cmap=pltsettings[0])
            varlabel = pltsettings[1]
            varrange = pltsettings[2]

            if projection == '2dxz':
                axs[t,v].set_yticks(y)
                axs[t,v].set_ylabel('z [m agl]')
            #plt.clim(varrange[0],varrange[1])
            #plt.axis('off')
            #ax.set_xlabel('x [m]')
        fig.colorbar(im, ax=axs[t,v], orientation='horizontal', label = varlabel)
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.{projection}-{idx}.png', dpi = 300)
    plt.close()


################################################################""

pfb_outputvariables = ['press','satur'] #,'velx','vely','velz']
plot_projections = ['2dxy','2dxz','3d']

nt = list_dumptimes()
varname = 'satur'
var_4D = readpfblist_to_4Darray(path,runname,varname)
#print(press_4D.shape)

yidx = 0
zidx = data.shape[0]-1
xidx = int(data.shape[2]/2)

#plot2dxz(var_4D,yidx,varname)
#plot2dxy(var_4D,zidx-1,varname)

plotmosaic(runname,'2dxz',yidx)
plotmosaic(runname,'2dxy',zidx)
#plotmosaic(runname,'2dyz',xidx)
