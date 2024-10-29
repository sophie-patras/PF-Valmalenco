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
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'

###############################################################################################
# INPUT

path = '/home/patras/PF-Valmalenco/outputs/'
foldername = 'CLM_V52'
runname = 'CLM_V5'

###############################################################################################

run = Run.from_definition(f'{path}{foldername}/{runname}.pfidb')
data = run.data_accessor

dx = data.dx
dy = data.dy
dz = data.dz
print(dz)

nx = data.shape[2]
ny = data.shape[1]
nz = data.shape[0]
#nt = data.time

porosity = data.computed_porosity
specific_storage = data.specific_storage
mask = data.mask
slopex = data.slope_x               # shape (ny, nx)
slopey = data.slope_y               # shape (ny, nx)
mannings = data.mannings

def read_pfidb(path,runname):
    
    run = Run.from_definition(f'{path}{foldername}/{runname}.pfidb')
    data = run.data_accessor

    dx = data.dx
    dy = data.dy
    dz = data.dz
    print(dz)

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
    print(centereddepth)
    return centereddepth

def list_dumptimes(): #(path,runname):
    ndt = 7
    return ndt

def readpfblist_to_4Darray(path,runname):

    nt = list_dumptimes() #path,runname)
    variable4D_rawpfb = np.empty((nt,nz,ny,nx))
    for t in range(nt):
        variable4D_rawpfb[t] = data.pressure
        data.time +=1

    return variable4D_rawpfb

def plot2dxz(variable4D_rawpfb,ny):

    array = variable4D_rawpfb[:,:,ny,:]
    nt = array.shape[0]

    x = np.linspace(0,nx-1,nx)
    z = layers_centereddepth()
    z = [round(z[i],4) for i in range(nz)]

    for t in range(nt):

        fig, ax = plt.subplots()
        plt.pcolormesh(x, z, array[t], cmap='Blues')
        varlabel = 'h [m above layer]'
        varrange = [array.min(), array.max()]

        #plt.imshow(data,cmap = colorofmap, origin="lower", extent=[0,10,-1,0], aspect=2) #, interpolation='nearest'/'none')
        
        ax.set(yticks=z, yticklabels=z)
        plt.colorbar(label = varlabel)
        #plt.colorbar(orientation='vertical', label = r'pressure head vs sl1 [m]', fraction = 0.04, pad = 0.0)
        plt.clim(varrange[0],varrange[1])
        
        #plt.axis('off')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('z [m agl]')
        plt.savefig(f'{path_fig}{foldername}.press.{str(t).zfill(5)}.2dxz.png', dpi = 300)
        plt.close()

def plot2dxy(variable4D_rawpfb,nz):

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
        plt.savefig(f'{path_fig}{foldername}.press.{str(t).zfill(5)}.2dxy-test.png', dpi = 300)
        plt.close()

def variablescaling(variable):

    if variable[:-1]=='vel':
        colorofmap = 'Spectral'
        varlabel = '$v_' + variable[-1] + '$ [m/h]'
        varrange = [data.min(),data.max()]
    elif variable == 'vel':
        colorofmap = 'Spectral'
        varlabel = 'v [m/h]'
        varrange = [data.min(),data.max()]
    elif variable == 'press':
        colorofmap = 'Blues'
        varlabel = 'h [m above layer]'
        varrange = [data.min(), data.max()] #min(data.max(),2)]
        #varrange = [-0.7, 0.1]
    elif variable == 'satur':
        colorofmap = 'Blues'
        varlabel = 'S [-]'
        #varrange = [data.min(), data.max()]
        varrange = [0, 1]

    return colorofmap, varlabel, varrange


press_4D = readpfblist_to_4Darray(path,runname)
print(press_4D.shape)
plot2dxz(press_4D,0)
plot2dxy(press_4D,1)

