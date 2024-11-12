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
from read_log import *

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant

###############################################################################################
# INPUT

path = '/home/patras/PF-Test/DumbSquare/outputs/'

foldername = 'DS.c100s1_v28' #'DS.c100s1_v24'
runname = 'DS.c100s1'

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

dt_real = dump_to_simulatedtimes_equivalent(path,foldername,runname)
print(f'dumptimes {dt_real}')
nt = len(dt_real)

porosity = data.computed_porosity
specific_storage = data.specific_storage
mask = data.mask
slopex = data.slope_x               # shape (ny, nx)
slopey = data.slope_y               # shape (ny, nx)
mannings = data.mannings

press_BC = data.pressure_boundary_conditions
print(press_BC)
# data.flow_boundary_conditions doesn't exist.


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

#def list_dumptimes(): #(path,runname):
#    ndt = 7 + 1
#    return ndt

#################################
# VARIABLES

def readpfblist_to_3Dtarray(path,runname,variable):
    # read raw pfb output 3D data 3D for each timestep, to 4D array
    data.time=0
    if variable == 'press':
        variable3Dt = np.empty((nt,nz,ny,nx))
        for t in range(nt):
            variable3Dt[t] = data.pressure
            data.time +=1
    if variable == 'satur':
        variable3Dt = np.empty((nt,nz,ny,nx))
        for t in range(nt):
            variable3Dt[t] = data.saturation
            data.time +=1
    if variable == 'velx':
        variable3Dt = np.empty((nt,nz,ny,nx+1))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable3Dt[t] = read_pfb(filename)
    if variable == 'vely':
        variable3Dt = np.empty((nt,nz,ny+1,nx))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable3Dt[t] = read_pfb(filename)
    if variable == 'velz':
        variable3Dt = np.empty((nt,nz+1,ny,nx))
        for t in range(nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable3Dt[t] = read_pfb(filename)
    if variable == 'subsurface_storage':
        variable3Dt = np.empty((nt,nz,ny,nx))
        for t in range(nt):
            variable3Dt[t] = data.subsurface_storage
            data.time +=1
    if variable == 'et':
        variable3Dt = np.empty((nt,nz,ny,nx))
        for t in range(nt):
            variable3Dt[t] = data.et
            data.time +=1

    return variable3Dt

def readpfblist_to_2Dtarray(path,runname,variable):
    # elaborate 2D(t) arrays, from pressure, saturation raw results
    if variable == 'surface_storage':
        variable2Dt = np.empty((nt,ny,nx))
        for t in range(nt):
            variable2Dt[t] = data.surface_storage
            data.time +=1
    if (variable == 'overlandsum') | (variable == 'overland_bc_flux'):
        variable2Dt = np.zeros((nt,ny,nx))
        for t in range(1,nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable2Dt[t] = read_pfb(filename)
    if variable == 'water_table':
        variable2Dt = np.zeros((nt,ny,nx))
        # or data.wtd
        saturation_3Dt = readpfblist_to_3Dtarray(path,runname,'satur')
        pressure_3Dt = readpfblist_to_3Dtarray(path,runname,'press')
        for t in range(1,nt):
            variable2Dt[t] = calculate_water_table_depth(pressure_3Dt[t], saturation_3Dt[t], dz)
    if variable == 'overland_flow':
        # flow_method='OverlandKinematic', default
        variable2Dt = np.zeros((nt,ny,nx))
        pressure_3Dt = readpfblist_to_3Dtarray(path,runname,'press')
        for t in range(1,nt):
            variable2Dt[t] = calculate_overland_flow_grid(pressure_3Dt[t], slopex, slopey, mannings, dx, dy, epsilon=1e-5, mask=mask)

    return variable2Dt

def discharge_3Dt(direction):
    # velocity to discharge (3D arrays), along 1 chosen direction
    if direction == 'y':
        array_3Dt = readpfblist_to_3Dtarray(path,runname,'vely')
        for z in range(nz):
            array_3Dt[:,z] = array_3Dt[:,z]*dx*dz[z]
    if direction == 'z':
        array_3Dt = readpfblist_to_3Dtarray(path,runname,'velz')
        array_3Dt = array_3Dt *dx*dy
    if direction == 'x':
        array_3Dt = readpfblist_to_3Dtarray(path,runname,'velx')
        for z in range(nz):
            array_3Dt[:,z] = array_3Dt[:,z]*dy*dz[z]

    return array_3Dt

def layer_mean1Dt(array_2Dt):
    # average a 2D array for each timestep
    array_1Dt = np.zeros(nt)
    for t in range(nt):
        array_1Dt[t] = np.mean(array_2Dt[t])
    return array_1Dt

def layer_sum1Dt(array_2Dt):
    # average a 2D array for each timestep
    array_1Dt = np.zeros(nt)
    for t in range(nt):
        array_1Dt[t] = np.sum(array_2Dt[t])
    return array_1Dt

#################################################################################"
# PLOTTING SETTINGS

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
    elif (variable == 'overlandsum') | (variable == 'overland_bc_flux'):
        colorofmap = 'jet'
        varlabel = 'overland flux [m$^3$/h]'
        varrange = [array.min(), array.max()] #array.min(),array.max()]
    elif variable == 'water_table':
        colorofmap = 'Blues'
        varlabel = 'water table [m bgl]'
        varrange = [array.min(), array.max()]
    elif variable == 'overland_flow':
        colorofmap = 'jet'
        varlabel = 'overland flow [m$^3$/h]'
        varrange = [array.min(), array.max()]

    return colorofmap, varlabel, varrange

################################################################################
# PLOT PAINTING/MOSAIC

def plotsingle_proj2d(variable4D_rawpfb,varname,dt,projection,idx):

    array = variable4D_rawpfb[:,:,ny,:]

    x = np.linspace(0,nx-1,nx)
    z = z_centers

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

def plotmosaic_proj2d(runname,projection,idx):

    pfb_outvariables = pfb_out3dtvariables
    if (projection == '2dxy') & (idx == nz-1):
        pfb_outvariables = pfb_out3dtvariables + pfb_out2dtvariables
    
    nbvar = len(pfb_outvariables)

    fig, axs = plt.subplots(nt,nbvar,figsize=(nbvar*6,nt*5),
                            gridspec_kw={'width_ratios': [1]*nbvar,
                                        'height_ratios': [1]*nt}) #
                                        #'wspace': 0.4,
                                        #'hspace': 0.4})

    for v in range(nbvar):
        vname = pfb_outvariables[v]
        print(vname, ' axs loop')

        projectioninfo = projscaling(vname, projection)
        x = projectioninfo[0]
        y = projectioninfo[1]
        xtit = projectioninfo[2]
        ytit = projectioninfo[3]
        xlim = projectioninfo[4]
        ylim = projectioninfo[5]

        if vname in pfb_out3dtvariables:
            array_3Dt = readpfblist_to_3Dtarray(path,runname,vname) #3Dt, future 2Dt
            array_2Dt = projected_array(array_3Dt,projection,idx)
        elif vname in pfb_out2dtvariables:
            array_2Dt = readpfblist_to_2Dtarray(path,runname,vname) #3Dt, future 2Dt
        #print(array_2Dt.shape)
        
        pltsettings = variablescaling(array_2Dt, vname)
        varcolor = pltsettings[0]
        varlabel = pltsettings[1]
        varrange = pltsettings[2]

        #array_fantom = np.empty(array_2Dt[0].shape)

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
            realt = dt_real[t]
            axs[t, 0].set_title(f't={realt}h',loc='left') #, ha='center', va='center', transform=axs[t, v].transAxes)
            fig.colorbar(im, ax=axs[t,v], orientation='vertical') #, fraction=0.5, pad=0.04)

        axs[0,v].set_title(f'{varlabel}')
        
        #im_fantom = axs[t+1,v].pcolormesh(x,y,array_fantom)
        #axs[t+1,v].set_aspect('0.001')
        #axs[t+1,v].axis('off')
        #fig.colorbar(im, ax=axs[t+1,v], orientation='horizontal', label = varlabel, fraction=0.5, pad=0.04)
    
    #plt.title(f'Sequence of PF direct outputs for run {foldername}, in plan {projection}')
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.{projection}{idx}.png', dpi = 300)
    plt.close()

def plot_multiXt(runname, Xcp):

    X = Xcp[0]
    nbX = X.shape[0]
    loc_cp = Xcp[1]
    position_layer = nz-1

    tarray = dt_real #np.linspace(0,nt-1,nt)
    
    nbvar = len(pfb_out3dtvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    for v in range(nbvar):
        vname = pfb_out3dtvariables[v]
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        varlabel = pltsettings[1]

        for x in range(nbX):
            array_Xt = array_3Dt[:,position_layer,X[x,0],X[x,1]]
            axs[v].plot(tarray,array_Xt,label=loc_cp[x]+f'({position_layer},{X[x,0]},{X[x,1]})', marker='.')
    
        axs[v].grid(True)
        axs[v].legend()
        axs[v].set_xlabel('t [h]')
        axs[v].set_ylabel(f'{varlabel}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.X.png', dpi = 300)

def plot_zXt(runname, Xcp):

    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    tarray = dt_real #np.linspace(0,nt-1,nt)
    
    nbvar = len(pfb_out3dtvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    for v in range(nbvar):
        vname = pfb_out3dtvariables[v]
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        varlabel = pltsettings[1]

        for z in range(nz):
            array_Xt = array_3Dt[:,z,X[0],X[1]]
            axs[v].plot(tarray,array_Xt,label=f'layer {z}', marker='.')
    
        axs[v].grid(True)
        axs[v].legend()
        axs[v].set_xlabel('t [h]')
        axs[v].set_ylabel(f'{varlabel}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.varall.dtall.{Xcp[1]}z.png', dpi = 300)

####################################
# CONVERGENCE PROXY

def plot_compared_press(runname, Xcp):
    # compare pressure-assumed linear, vs, water table, (vs BC if point on y-lower)
    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    tarray = dt_real #np.linspace(0,nt-1,nt)
    
    nbvar = 1
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    # pressure head to press
    vname = 'press'
    v = 1
    array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    pltsettings = variablescaling(array_3Dt, vname)
    #varlabel = pltsettings[1]

    for z in range(nz):
        array_Xt = array_3Dt[:,z,X[0],X[1]] + z_centers[z]
        axs.plot(tarray,array_Xt,label=f'press(layer={z})', marker='.')

    # Water table
    vname = 'water_table'
    v = 2
    array_3Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
    pltsettings = variablescaling(array_3Dt, vname)
    #varlabel = pltsettings[1]
    array_Xt = array_3Dt[:,X[0],X[1]]
    axs.plot(tarray,array_Xt,label='pf-wt', marker='.')
    
    """
    # PBC y-lower
    v = 3
    ylower = press_BC['y-lower__alltime'] * np.ones(dt_real.shape)
    axs.plot(tarray,ylower,label='PBC y-lower')
    """

    axs.grid(True)
    axs.legend()
    axs.set_xlabel('t [h]')
    axs.set_ylabel('pressure [m agl]')

    plt.title('Comparison of linear pressure, water table (-1), press BC')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.presscomp.dtall.{loc_cp}.png', dpi = 300)

def plot_flow_at_boundaries():
    # for all faces in press_BC list - for CV sum should go to 0 ('balanced model')
    # normal to face is directed outward of domain

    sum = np.zeros(nt)

    fig, ax = plt.subplots(1,1,figsize=(6,5))
    # 'x-lower'
    array_Xt = (-1)*layer_sum1Dt(projected_array(discharge_3Dt('x'),'2dyz',0))
    sum = array_Xt
    ax.plot(dt_real,array_Xt,label='x-lower', marker='.')
    # 'x-upper'
    array_Xt = layer_sum1Dt(projected_array(discharge_3Dt('x'),'2dyz',nx))
    sum += array_Xt
    ax.plot(dt_real,array_Xt,label='x-upper', marker='.')
    # 'y-lower'
    array_Xt = (-1)*layer_sum1Dt(projected_array(discharge_3Dt('y'),'2dxz',0))
    sum += array_Xt
    ax.plot(dt_real,array_Xt,label='y-lower', marker='.')
    # 'y-upper'
    array_Xt = layer_sum1Dt(projected_array(discharge_3Dt('y'),'2dxz',ny))
    sum += array_Xt
    ax.plot(dt_real,array_Xt,label='y-upper', marker='.')
    # 'z-lower'
    array_Xt = (-1)*layer_sum1Dt(projected_array(discharge_3Dt('z'),'2dxy',0))
    sum += array_Xt
    ax.plot(dt_real,array_Xt,label='z-lower', marker='.')
    # 'z-upper'
    array_Xt = layer_sum1Dt(projected_array(discharge_3Dt('z'),'2dxy',nz))
    sum += array_Xt
    ax.plot(dt_real, array_Xt, label='z-upper', marker='.')
    # 'sum'
    #ax.plot(dt_real,sum,label='sum', color='k',marker='.')

    ax.grid(True)
    ax.legend()
    ax.set_xlabel('t [h]')
    ax.set_ylabel('boundary outflow [m$^3$/h]')

    plt.title('Total outflow at each boundary face and sum')
    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.outflow.Bf.dtall.png', dpi = 300)


################################################################""

pfb_out3dtvariables = ['press','satur','velx','vely','velz'] #'subsurface_storage'
pfb_out2dtvariables = ['water_table', 'overland_flow'] #'surface_storage', 'surface_storage', 'overland_bc_flux', 'overlandsum',
faces = ['x-lower', 'x-upper', 'y-lower', 'y-upper', 'z-lower', 'z-upper']
plot_projections = ['2dxy','2dxz','2dyz','1d','3d']

#varname = 'satur'
#var_4D = readpfblist_to_4Darray(path,runname,varname)
#print(press_4D.shape)

zidx = nz-1
yidx = 0
xidx = int(nx/2)

#plot2dxz(var_4D,yidx,varname)
#plot2dxy(var_4D,zidx-1,varname)

#plotmosaic_proj2d(runname,'2dxz',yidx)
#plotmosaic_proj2d(runname,'2dxy',zidx)
#plotmosaic_proj2d(runname,'2dyz',xidx)

#f_cp = '/home/patras/PF-Valmalenco/data/controlpoints.txt'
#Xidx_cp = read_cpcsv(f_cp)
#print(Xidx_cp)
#XP = [Xidx_cp[0][1],[Xidx_cp[1][1]]]

Xidx_cp = [np.array([[0,5],[5,5],[9,2],[9,5]]),['P1','P2','P3','P4']]
XP_center = [np.array([5,5]),['P2']]
XP_ylower = [np.array([0,5]),['P1']]
#plot_multiXt(runname,Xidx_cp)
#plot_zXt(runname,XP_center)
#plot_compared_press(runname,XP_ylower)

plot_flow_at_boundaries()