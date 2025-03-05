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
import matplotlib.colors as colors
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

###
plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#np.set_printoptions(precision=6)
###

#parser = argparse.ArgumentParser()
#parser.add_argument('filename', type = str)
#args = parser.parse_args()

#### INPUT FILES
#filename = args.filename
# SOILGRIDS
#path_in = '/home/patras/PF-Valmalenco/data/pfdata/' #  from input file folder
#filename = 'KS.c250.v8-331.pfb'

# COSMO
path_in = '/home/patras/PF-Valmalenco/data/pfdata/clm.v12-331/'
#filename = 'clm.v12-331.APCP.000121_to_000144.pfb'
vname = 'APCP'
filename = f'clm.v12-331.{vname}.000001_to_000024.pfb'

f_in = path_in + filename

sg_variables = ['KS','MRC_ths','MRC_thr','MRC_alp','MRC_n']
clm_variables = ['DSWR', 'DLWR', 'Press', 'SPFH', 'Temp',  'APCP', 'UGRD', 'VGRD']


####
data = read_pfb(f_in)
datashape = data.shape
#data = np.flip(data,axis=1)
print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)

# Apply mask where data <-3.1e+38 (ie. nodata_value)
#data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)

dx = 250.0
dy = 250.0
dz = [24.0, 12.0, 6.0, 3.0, 1.5, 0.7, 0.35, 0.25, 0.0175, 0.025, 0.025]

nx = data.shape[2]
ny = data.shape[1]
nz = data.shape[0]

x_centers = dx*(np.linspace(1,nx,nx)-1/2)
y_centers = dy*(np.linspace(1,ny,ny)-1/2)
x_faces = dx*(np.linspace(0,nx,nx+1))
y_faces = dy*(np.linspace(0,ny,ny+1))

#### PLOT SETTINGS

def projscaling(varname, projection):

    if projection == '2dxz':
        xlabels = 'x [m]'
        ylabels = 'z [m agl]'
        xlim = x_faces
        ylim = z_faces
    if projection == '2dxy':
        xlabels = 'x [m]'
        ylabels = 'y [m]'
        xlim = x_faces
        ylim = y_faces
    if projection == '2dyz':
        xlabels = 'y [m]'
        ylabels = 'z [m agl]'
        xlim = y_faces
        ylim = z_faces

    # velocities//i : faces centered values//i
    if varname == 'velx':
        if projection == '2dxz':
            x = x_faces
            y = z_centers
        if projection == '2dxy':
            x = x_faces
            y = y_centers
        if projection == '2dyz':
            x = y_centers
            y = z_centers

    if varname == 'vely':
        if projection == '2dxz':
            x = x_centers
            y = z_centers
        if projection == '2dxy':
            x = x_centers
            y = y_faces
        if projection == '2dyz':
            x = y_faces
            y = z_centers

    if varname == 'velz':
        if projection == '2dxz':
            x = x_centers
            y = z_faces
        if projection == '2dxy':
            x = x_centers
            y = y_centers
        if projection == '2dyz':
            x = y_centers
            y = z_faces

    if varname not in ['velx','vely','velz']:
        #print(f'{varname}')
    # cell centered values
        if projection == '2dxz':
            #print('test :',varname)
            x = x_centers
            y = z_centers
        if projection == '2dxy':
            x = x_centers
            y = y_centers
        if projection == '2dyz':
            x = y_centers
            y = z_centers

    return x,y,xlabels,ylabels,xlim,ylim

def variablescaling(array,variable):

    norm = None

    if variable=='APCP':
        colorofmap = 'Blues'
        varlabel = 'rainfall [mm/h]'
        varrange = [array.min(),array.max()]
        #norm = colors.LogNorm(vmin=array.min(), vmax=array.max())
    elif variable=='SPFH':
        colorofmap = 'jet'
        varlabel = r'SPFH [-]'
        varrange = [array.min(),array.max()]
    elif variable=='UGRD':
        colorofmap = 'jet'
        varlabel = r'$v_{x,wind}$ [m/h]'
        varrange = [array.min(),array.max()]
    elif variable=='VGRD':
        colorofmap = 'jet'
        varlabel = r'$v_{y,wind}$ [m/h]'
        varrange = [array.min(),array.max()]
    elif variable=='Temp':
        colorofmap = 'jet'
        varlabel = r'T [K]'
        varrange = [array.min(),array.max()]
    elif variable=='DSWR':
        colorofmap = 'jet'
        varlabel = r'DSWR [W/m$^2$]'
        varrange = [array.min(),array.max()]
    elif variable=='DLWR':
        colorofmap = 'jet'
        varlabel = r'DLWR [W/m$^2$]'
        varrange = [array.min(),array.max()]
    elif variable=='Press':
        colorofmap = 'jet'
        varlabel = r'atm press [Pa]'
        varrange = [array.min(),array.max()]

    elif variable == 'slopex':
        colorofmap = 'PuOr'
        varlabel = r"$i_x$ [-]"
        varrange = [array.min(),array.max()]
    elif variable == 'slopey':
        colorofmap = 'PuOr'
        varlabel = r"$i_y$ [-]"
        varrange = [array.min(),array.max()]
    elif variable == 'DEM':
        colorofmap = 'terrain'
        varlabel = r"$z$ [m asl]"
        varrange = [0, 4000] #[array.min(),array.max()] # [0, 4050]
    elif variable == 'KS':
        colorofmap = 'jet'
        varlabel = r'$K_{sat}$ [m/h]'
        varrange = [array.min(), array.max()]
        norm = colors.LogNorm(vmin=array.min(), vmax=array.max())
    elif variable == 'alpha':
        colorofmap = 'viridis'
        varlabel = r'$\alpha$ [m$^{-1}$]'
        varrange = [array.min(),array.max()]      
    elif variable == 'n':
        colorofmap = 'viridis'
        varlabel = r'$n$ [-]'
        varrange = [array.min(),array.max()] 
    elif variable == 'ssat':
        colorofmap = 'Blues'
        varlabel = r'$S_{sat}$ [-]'
        varrange = [array.min(),array.max()] 
    elif variable == 'sres':
        colorofmap = 'Blues'
        varlabel = r'$S_{res}$ [-]'
        varrange = [array.min(),array.max()]   
    else: # default
        colorofmap = 'Spectral_r'
        varlabel = f'{variable}'
        varrange = [array.min(),array.max()]

    return colorofmap, varlabel, varrange, norm

pltsettings = variablescaling(data, vname)
colorofmap = pltsettings[0]
varlabel = pltsettings[1]
varrange = pltsettings[2]
norm = pltsettings[3]

projectioninfo = projscaling(vname, '2dxy')
x = projectioninfo[0]
y = projectioninfo[1]
xtit = projectioninfo[2]
ytit = projectioninfo[3]

fig, axs = plt.subplots(datashape[0],1,figsize=(6,datashape[0]*5))

for idx in range(datashape[0]):
    print('layer',idx)
    array_2Dl = data[idx]
    #data = data_masked[p]
    #print(datap)
    if vname in sg_variables:
        layername = 'sl' + str(datashape[0]-idx)
    if vname in clm_variables:
        layername = f't={idx}h'
    #vmin = datap.min()
    #vmax = datap.max()
    #vmean = np.mean(datap)
    #print(vmin,vmean,vmax)

    #plt.imshow(datap,cmap = colorofmap, norm=norm) #, interpolation='nearest'
    # plt.colorbar(label = varlabel)
    im = axs[idx].pcolormesh(x, y, array_2Dl, shading='auto', cmap=colorofmap, vmin = varrange[0], vmax=varrange[1])
    fig.colorbar(im, ax=axs[idx], orientation='vertical') #, fraction=0.5, pad=0.04)
    axs[idx].set_xlabel(xtit)
    axs[idx].set_ylabel(ytit)
    axs[idx].set_title(layername,loc='left') #, ha='center', va='center', transform=axs[t, v].transAxes)
    axs[idx].set_title(varlabel) #, ha='center', va='center', transform=axs[t, v].transAxes)
    #plt.show()

plt.savefig('/mnt/c/Users/Sophie/Documents/4-Figures/'+filename+'.2dxy.png')
plt.close()
