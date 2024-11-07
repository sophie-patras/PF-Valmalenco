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

from Xcrs_to_Xidx import *

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant

###############################################################################################
# INPUT

path = '/home/patras/PF-Test/DumbSquare/outputs/'
foldername = 'DS.c1_v7R'
runname = 'DS.c1'

"""
path = '/home/patras/PF-Valmalenco/outputs/'
foldername = 'CLM_V22'
runname = 'CLM_V2'
"""

#path = '/home/patras/PF-Test/LW/outputs/'
#foldername = 'LW_var_dz_spinup'
#runname = 'LW_var_dz_spinup'

###############################################################################################

run = Run.from_definition(f'{path}{foldername}/{runname}.pfidb')
data = run.data_accessor

#print(data)

dx = data.dx
dy = data.dy
dz = data.dz
#print(dz)

nx = data.shape[2]
ny = data.shape[1]
nz = data.shape[0]
#nt = data.time

# domain rect cells
x_centers = dx*(np.linspace(1,nx,nx)-1/2)
y_centers = dy*(np.linspace(1,ny,ny)-1/2)
x_faces = dx*(np.linspace(0,nx,nx+1))
y_faces = dy*(np.linspace(0,ny,ny+1))

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

#def layers_centereddepth():
    
dzscaled = np.concatenate((dz,[0]))
z_faces = (-1)*np.array([sum(dzscaled[i:]) for i in range(nz+1)])
print('z_faces',z_faces)
z_centers = (z_faces[:-1]+z_faces[1:])/2
print('z_centers',z_centers)
z_centers = [round(z_centers[i],4) for i in range(nz)]
#return z_centers, z_faces

def list_dumptimes(): #(path,runname):
    ndt = 7 + 1
    return ndt

def readpfblist_to_3Dtarray(path,runname,variable):

    nt = list_dumptimes() #path,runname)
    data.time=0
    if variable == 'press':
        variable4D_rawpfb = np.empty((nt,nz,ny,nx))
        for t in range(nt):
            variable4D_rawpfb[t] = data.pressure
            data.time +=1
    if variable == 'satur':
        variable4D_rawpfb = np.empty((nt,nz,ny,nx))
        for t in range(nt):
            variable4D_rawpfb[t] = data.saturation
            data.time +=1
    if variable == 'velx':
        variable4D_rawpfb = np.empty((nt,nz,ny,nx+1))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable4D_rawpfb[t] = read_pfb(filename)
    if variable == 'vely':
        variable4D_rawpfb = np.empty((nt,nz,ny+1,nx))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable4D_rawpfb[t] = read_pfb(filename)
    if variable == 'velz':
        variable4D_rawpfb = np.empty((nt,nz+1,ny,nx))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable4D_rawpfb[t] = read_pfb(filename)
    if variable == 'overlandsum':
        variable4D_rawpfb = np.empty((nt,nz+1,ny,nx))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable4D_rawpfb[t] = read_pfb(filename)
    if variable == 'overland_bc_flux':
        variable4D_rawpfb = np.empty((nt,nz+1,ny,nx))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable4D_rawpfb[t] = read_pfb(filename)

    return variable4D_rawpfb


def projected_array(array_3Dt,projection,idx):

    array = array_3Dt
    if projection == '2dxz':
        array_2Dt = array[:,:,idx]
    if projection == '2dxy':
        array_2Dt = array[:,idx]
    if projection == '2dyz':
        array_2Dt = array[:,:,:,idx]

    return array_2Dt

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

    if (varname == 'press') | (varname == 'satur'):
    # cell centered values
        if projection == '2dxz':
            print('test :',varname)
            x = x_centers
            y = z_centers
        if projection == '2dxy':
            x = x_centers
            y = y_centers
        if projection == '2dyz':
            x = y_centers
            y = z_centers

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

    return x,y,xlabels,ylabels,xlim,ylim

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
        varrange = [array.min(), array.max()] #min(array.max(),2)]
        #varrange = [-0.7, 0.1]
    elif variable == 'satur':
        colorofmap = 'Blues'
        varlabel = 'S [-]'
        #varrange = [data.min(), data.max()]
        varrange = [max(0, array.min()), min(array.max(), 1)]

    return colorofmap, varlabel, varrange

def plotsingle_proj2d(variable4D_rawpfb,varname,dt,projection,idx):

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

def plotmosaic_proj2d(runname,projection,idx):

    nbvar = len(pfb_outputvariables)

    fig, axs = plt.subplots(nt+1,nbvar,figsize=(nbvar*6,nt*5),
                            gridspec_kw={'width_ratios': [1, 1, 1, 1, 1],
                                        'height_ratios': [1, 1, 1, 1, 1, 1, 1, 1, 0.5]}) #
                                        #'wspace': 0.4,
                                        #'hspace': 0.4})

    for v in range(nbvar):
        vname = pfb_outputvariables[v]
        print(vname, ' axs loop')

        projectioninfo = projscaling(vname, projection)
        x = projectioninfo[0]
        y = projectioninfo[1]
        xtit = projectioninfo[2]
        ytit = projectioninfo[3]
        xlim = projectioninfo[4]
        ylim = projectioninfo[5]

        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        array_2Dt = projected_array(array_3Dt,projection,idx)
        pltsettings = variablescaling(array_2Dt, vname)
        varcolor = pltsettings[0]
        varlabel = pltsettings[1]
        varrange = pltsettings[2]

        array_fantom = np.empty(array_2Dt[0].shape)

        for t in range(nt):
            print(t, 't axs loop', np.mean(array_2Dt[t]))
            im = axs[t,v].pcolormesh(x, y, array_2Dt[t], shading='auto', cmap=varcolor, vmin = varrange[0], vmax=varrange[1] )

            if (projection == '2dxz') | (projection == '2dyz'):
                axs[t,v].set_yticks(y)
                #axs[t,v].set_aspect(10000)
            elif projection == '2dxy':
                axs[t,v].set_aspect('equal')
            
            axs[t,v].set_xlabel(xtit)
            #axs[t,v].set_xlim(xlim[0],xlim[-1])
            axs[t,v].set_ylabel(ytit)
            #axs[t,v].set_ylim(ylim[0],ylim[-1])

            #plt.axis('off')

            axs[t, 0].set_title(f't={t}d',loc='left') #, ha='center', va='center', transform=axs[t, v].transAxes)
        axs[0,v].set_title(f'{varlabel}')
        #axs[0,0].set_title(f't=0d \t\t\t {varlabel}',loc='left')
        im_fantom = axs[t+1,v].pcolormesh(x,y,array_fantom)
        axs[t+1,v].set_aspect('0.001')
        axs[t+1,v].axis('off')
        fig.colorbar(im, ax=axs[t+1,v], orientation='horizontal', label = varlabel, fraction=0.5, pad=0.04)
        
    #plt.title(f'Sequence of PF direct outputs for run {foldername}, in plan {projection}')
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.{projection}-{idx}.png', dpi = 300)
    plt.close()

def plot_multiXt(runname, Xcp):

    X = Xcp[0]
    nbX = X.shape[0]
    loc_cp = Xcp[1]
    position_layer = nz-1
    tarray = np.linspace(0,nt-1,nt)
    
    nbvar = len(pfb_outputvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    for v in range(nbvar):
        vname = pfb_outputvariables[v]
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        varlabel = pltsettings[1]

        for x in range(nbX):
            array_Xt = array_3Dt[:,position_layer,X[x,0],X[x,1]]
            axs[v].plot(tarray,array_Xt,label=loc_cp[x]+f'({position_layer},{X[x,0]},{X[x,1]})')
    
        axs[v].grid(True)
        axs[v].legend()
        axs[v].set_xlabel('t [d]')
        axs[v].set_ylabel(f'{varlabel}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.X.png', dpi = 300)

def plot_zXt(runname, Xcp):

    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    tarray = np.linspace(0,nt-1,nt)
    
    nbvar = len(pfb_outputvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    for v in range(nbvar):
        vname = pfb_outputvariables[v]
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        varlabel = pltsettings[1]

        for z in range(nz):
            array_Xt = array_3Dt[:,z,X[0],X[1]]
            axs[v].plot(tarray,array_Xt,label=f'layer {z}')
    
        axs[v].grid(True)
        axs[v].legend()
        axs[v].set_xlabel('t [d]')
        axs[v].set_ylabel(f'{varlabel}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.{Xcp[1]}z.png', dpi = 300)


################################################################""

pfb_outputvariables = ['press','satur','velx','vely','velz'] #,'overland_bc_flux', 'overlandsum']
plot_projections = ['2dxy','2dxz','2dyz','1d','3d']

nt = list_dumptimes()
#varname = 'satur'
#var_4D = readpfblist_to_4Darray(path,runname,varname)
#print(press_4D.shape)

zidx = data.shape[0]-1
yidx = 0
xidx = int(data.shape[2]/2)

#plot2dxz(var_4D,yidx,varname)
#plot2dxy(var_4D,zidx-1,varname)

plotmosaic_proj2d(runname,'2dxz',yidx)
plotmosaic_proj2d(runname,'2dxy',zidx)
plotmosaic_proj2d(runname,'2dyz',xidx)

#f_cp = '/home/patras/PF-Valmalenco/data/controlpoints.txt'
#Xidx_cp = read_cpcsv(f_cp)
#print(Xidx_cp)
#XP = [Xidx_cp[0][1],[Xidx_cp[1][1]]]

Xidx_cp = [np.array([[1,5],[5,5],[9,5],[9,2]]),['P1','P2','P3','P4']]
XP = [np.array([5,5]),['P2']]

plot_multiXt(runname,Xidx_cp)
plot_zXt(runname,XP)