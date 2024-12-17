#!/usr/bin/python3
# -*- coding: utf-8 -*-
#%%

"""
For graphical output
Read all times <runname>.out.<variable>.<dumptime>.pfb files (p-layers) of an assigned folder path
Functions : 
- 3d
- 2dproj (xy,xz,yz) @ given missing index
- 1dz(t)
...etc.

Created on Mon Oct 28 16:41:20 2024
@author: S.P.


version : DRAFT v1
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors # using LogNorm
from mpl_toolkits.mplot3d import Axes3D
# import imageio

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb
from parflow import Run
from parflow.tools.hydrology import calculate_surface_storage, calculate_water_table_depth, calculate_overland_flow_grid, compute_hydraulic_head, calculate_overland_fluxes, calculate_subsurface_storage

# check what is the difference btw compute_water_table_depth and calculate_water_table_depth
# check compute_hydraulic_head

from Xcrs_to_Xidx import *
from read_log import *


# coment if not defined
plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')

#####################################################################################"
# OUTPUT PATH

#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant

###############################################################################################
# INPUT PATH
# define the complete path to your output .pfidb and .pfb files
# such that <runname>.pfidb can be found in : f'{path}{foldername}{runname}.pfidb'

## DS
path = '/home/patras/PF-Test/DumbSquare/outputs/'
foldername = 'DS.c100s05_v35' #'DS.c100s1_v28'
runname = 'DS.c100s05'

## MX
#path = '/home/patras/PF-Test/Maxwell2013/outputs/'
#foldername = 'MX.c100s1bcx_v5'
#runname = 'MX.c1s1'

## VM
path = '/home/patras/PF-Valmalenco/outputs/'
foldername = 'CLM_V52'
runname = 'CLM_V5'

## parflow test cases
#path = '/home/patras/PF-Test/LW/outputs/'
#foldername = 'LW_var_dz_spinup'
#runname = 'LW_var_dz_spinup'

###############################################################################################

run = Run.from_definition(f'{path}{foldername}/{runname}.pfidb')
data = run.data_accessor

#print(data)

#########################################
# Discretization parameters
dx = data.dx
dy = data.dy
dz = data.dz
#print(dz)

nx = data.shape[2]
ny = data.shape[1]
nz = data.shape[0]

# matrix of dz
dz_3D = np.ones((nz,ny,nx))
for i in range(nz):
    dz_3D[i] = dz_3D[i]*dz[i]

#nt = data.time

X_origin = [0,0,0] # change with real CRS coordinate
print(X_origin[0])

# domain rect cells
x_centers = dx*(np.linspace(1,nx,nx)-1/2)
y_centers = dy*(np.linspace(1,ny,ny)-1/2)
x_faces = dx*(np.linspace(0,nx,nx+1))
y_faces = dy*(np.linspace(0,ny,ny+1))
#print(x_centers.shape)
#print(y_centers.shape)

dt_real = dump_to_simulatedtimes_equivalent(path,foldername,runname)
#print(f'dumptimes {dt_real}')
nt = len(dt_real)

#def layers_centereddepth():
dzscaled = np.concatenate((dz,[0]))
z_faces = (-1)*np.array([sum(dzscaled[i:]) for i in range(nz+1)])
print('z_faces',z_faces)
z_centers = (z_faces[:-1]+z_faces[1:])/2
print('z_centers',z_centers)
z_centers_rnd = [round(z_centers[i],4) for i in range(nz)]
print('z_centers_rnd',z_centers_rnd)

# matrix of z_agl
z_centers_3D = np.ones((nz,ny,nx))
for i in range(nz):
    z_centers_3D[i] = z_centers_3D[i]*z_centers[i]

mask = data.mask
slopex = data.slope_x               # shape (1,ny, nx)
slopey = data.slope_y               # shape (1,ny, nx)
#print(slopey[0,0,0])
#print(slopex.shape)

def compute_dem():
    array_2D = np.empty((ny,nx))
    array_2D[0,0] = X_origin[0] + slopex[0,0,0]*dx/2
    for i in range(1,nx):
        array_2D[0,i] = array_2D[0,i-1] + (slopex[0,0,i-1]+slopex[0,0,i])*dx/2
    for i in range(0,nx):
        for j in range(1,ny):
            array_2D[j,i] = array_2D[j-1,i] + (slopey[0,j-1,i]+slopey[0,j,i])*dy/2
    return array_2D

# surface dem [asl]
z_dem_2D = compute_dem()
#print(z_dem_2D)

# matrix of z_asl
def compute_zcentersasl_3D():
    array_3D =  z_centers_3D
    for k in range(nz):
        array_3D[k] += z_dem_2D
    return array_3D

z_centers_asl_3D = compute_zcentersasl_3D()

dt_real = dump_to_simulatedtimes_equivalent(path,foldername,runname)
#print(f'dumptimes {dt_real}')
nt = len(dt_real)

#################################
# DISTRIBUTED INPUT PARAMETERS

porosity = data.computed_porosity
specific_storage = data.specific_storage
mannings = data.mannings
permx = data.computed_permeability_x
permy = data.computed_permeability_y
permz = data.computed_permeability_z
perm_avg = np.mean(permz)
print(perm_avg)

press_BC = data.pressure_boundary_conditions
print('Pressure head boundary conditions : ', press_BC)
# data.flow_boundary_conditions doesn't exist.

#################################
# VARIABLES

def readpfblist_to_3Dtarray(path,runname,variable):
    # read raw pfb output 3D data 3D for each timestep, to 4D array
    data.time=0
    if variable == 'press':
        # pressure head (as gravity=1, density=1)
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
        for t in range(0,nt):
            filename = path + foldername + '/' + runname + ".out." + variable + "." + str(t).zfill(5) + ".pfb"
            variable2Dt[t] = read_pfb(filename)
    if variable == 'water_table':
        variable2Dt = np.zeros((nt,ny,nx))
        # or data.wtd
        saturation_3Dt = readpfblist_to_3Dtarray(path,runname,'satur')
        pressure_3Dt = readpfblist_to_3Dtarray(path,runname,'press')
        for t in range(0,nt):
            variable2Dt[t] = calculate_water_table_depth(pressure_3Dt[t], saturation_3Dt[t], dz)
    if variable == 'overland_flow':
        # flow_method='OverlandKinematic', default
        variable2Dt = np.zeros((nt,ny,nx))
        pressure_3Dt = readpfblist_to_3Dtarray(path,runname,'press')
        for t in range(0,nt):
            variable2Dt[t] = calculate_overland_flow_grid(pressure_3Dt[t], slopex, slopey, mannings, dx, dy, epsilon=1e-5, mask=mask)
        variable2Dt = variable2Dt #* 1/3600 # m^3/s

        #flow = calculate_overland_flow_grid(press, slopex, slopey, mannings, dx, dy, flow_method='OverlandKinematic', epsilon=1e-5, mask=mask)
        #flx, fly = calculate_overland_fluxes(press, slopex, slopey, mannings, dx, dy)

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
def dischargenorm_3Dt():
    # cell centered
    vel = velocitynorm_3Dt()
    q = vel * dz_3D * dy * dx
    return q

def velocitynorm_3Dt():
    
    velx_centered = face_to_centered_variable(readpfblist_to_3Dtarray(path, runname, 'velx'),'x')
    vely_centered = face_to_centered_variable(readpfblist_to_3Dtarray(path, runname, 'vely'),'y')
    velz_centered = face_to_centered_variable(readpfblist_to_3Dtarray(path, runname, 'velz'),'z')
    array_3Dt = vector_norm(velx_centered,vely_centered,velz_centered)
    return array_3Dt
  
def hydraulichead_3Dt():
    # centered variable
    # nb. Bernoulli assumptions is not verified (not Newtonian), thus h is an estimate of the hydraulic head
    vel = velocitynorm_3Dt()
    press = readpfblist_to_3Dtarray(path, runname, 'press') # pressure head
    # h = p + z + u^2/2
    # Element-wise operation
    h = press + z_centers_3D #+ vel**2
    return h

def hydraulichhead_3Dt_bis():
    press = readpfblist_to_3Dtarray(path, runname, 'press') # pressure head
    h = compute_hydraulic_head(press, 0, dz) # hh = hp + hz
    return h

#def overland_flow_manual():
    # Manning-strickler equation. V = KsRh^2/3i1/2
    # y direction

def centered_var_generator(vname):
# pfb_cellcentered3dtvariables = ['press', 'satur', 'velx_centers','vely_centers','velz_centers','vel_norm', 'hydraulic_head']

    if vname == 'v':
        array3Dt = velocitynorm_3Dt()
    elif vname == 'H':
        array3Dt = hydraulichead_3Dt()
    elif vname in pfb_out3dtvariables:
        array3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    elif vname in ['vx_c','vy_c','vz_c']:
        if vname == 'vx_c':
            array3Dt = readpfblist_to_3Dtarray(path,runname,'velx')
            array3Dt = face_to_centered_variable(array3Dt,'x')
        elif vname == 'vy_c':
            array3Dt = readpfblist_to_3Dtarray(path,runname,'vely')
            array3Dt = face_to_centered_variable(array3Dt,'y')
        elif vname == 'vz_c':
            array3Dt = readpfblist_to_3Dtarray(path,runname,'velz')
            array3Dt = face_to_centered_variable(array3Dt,'z')

    return array3Dt

#################################################################################
# Array operations

# Apply mask where data <-3.1e+38 (ie. nodata_value)
# data_masked = np.ma.masked_where(data_to_plot <= -3.e+38, data_to_plot, copy=True)

def face_to_centered_variable(array3Dt_faced,direction):
    # simplest linear interpolation from 2 closest point (average)
    # from face centered scheme to cell centered
    # upward scheme
    variable3Dt_centered = np.empty((nt,nz,ny,nx))
    if direction == 'x':
        variable3Dt_centered = array3Dt_faced[:,:,:,1:]-array3Dt_faced[:,:,:,:-1]
    if direction == 'y':
        variable3Dt_centered = array3Dt_faced[:,:,1:,:]-array3Dt_faced[:,:,:-1,:]
    if direction == 'z':
        variable3Dt_centered = array3Dt_faced[:,1:,:,:]-array3Dt_faced[:,:-1,:,:]

    return variable3Dt_centered

def vector_norm(datax, datay, dataz):
    #array_3Dt = np.sqrt(np.square(datax[:,:,:,:-1]) + np.square(datay[:,:,1:,:]) + np.square(dataz[:,:-1,:,:]))
    # all cell centered variables
    array_3Dt = np.sqrt(np.square(datax) + np.square(datay) + np.square(dataz))
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
# PLOT SETTINGS x variable

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

    norm = None

    if variable[:-1]=='vel':
        colorofmap = 'Spectral_r'
        varlabel = '$v_' + variable[-1] + '$ [m/h]'
        varrange = [array.min(),array.max()]
        #norm = colors.LogNorm(vmin=array.min(), vmax=array.max())
    elif variable == 'v':
        colorofmap = 'Spectral_r'
        varlabel = 'v [m/h]'
        varrange = [array.min(),array.max()]
        #norm = colors.LogNorm(vmin=array.min(), vmax=array.max())
    elif variable == 'press':
        colorofmap = 'Blues'
        varlabel = r'$p/\rho g$ [m above layer]'
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
        colorofmap = 'Spectral_r'
        varlabel = 'overland flow [m$^3$/h]'
        varrange = [array.min(), array.max()]
    elif variable == 'H':
        colorofmap = 'Blues'
        varlabel = 'hydraulic head [m agl.]'
        varrange = [array.min(),array.max()]
    else: # default
        colorofmap = 'Spectral_r'
        varlabel = f'{variable}'
        varrange = [array.min(),array.max()]

    return colorofmap, varlabel, varrange, norm

def tidxscaling(dt_fullreal,mod,**kwargs):
    # return idx of time
    nbofmosaic = kwargs.get('c',int)
    if mod == "all":
        idx = np.linspace(0,nt-1,nt)
        idx = idx.astype(int)
    if mod == "final":
        idx = [nt-1]
    elif mod == "intall":
        tinit = kwargs.get('d',float)
        tfin = kwargs.get('e',float)
        idxinit = np.where(dt_real<=tinit)[-1][-1]
        #print(idxinit[-1])
        idxfin = np.where(dt_real>=tfin)[-1][0]
        #print(idxfin[0])
        idx = np.linspace(idxinit,idxfin,idxfin-idxinit+1)
        #print(idx)
        idx = idx.astype(int)
    elif (mod!="all") and (mod!="intall") and (nbofmosaic != None) and (nt>(2*nbofmosaic-1)): # kwargs.get('c',float) = nb of dumptimes
        print('passed')
        if mod == "ed": # "equaldist_indumptimes"
            idx = np.linspace(0,nt-1,nbofmosaic)
            idx = idx.astype(int)
        if mod == "beg":
            idx = np.linspace(0,nbofmosaic-1,nbofmosaic)
            idx = idx.astype(int)
        if mod == "int": #interval of time
            tinit = kwargs.get('d',float)
            tfin = kwargs.get('e',float)
            idxinit = np.where(dt_real<=tinit)[-1][-1]
            idxfin = np.where(dt_real>=tfin)[-1][0]
            idx = np.linspace(idxinit,idxfin,nbofmosaic)
            idx = idx.astype(int)
        if mod == "fin":
            idx = np.linspace(nt-nbofmosaic,nt-1,nbofmosaic)
            idx = idx.astype(int)
    else:
        idx = np.linspace(0,nt-1,nt)
        idx = idx.astype(int)
    return idx

#dt_new = dtscaling(dt_real,"equaldist_indumptimes",5)
#print(dt_new)
################################################################################
# PLOT 'PAINTINGS'

# https://jackmckew.dev/3d-terrain-in-python.html

def plot_3d_geom():

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)
    
    #x,y = np.meshgrid(np.arange(data_3Dt.shape[3]), np.arange(data_3Dt.shape[2]))
    x_c, y_c = np.meshgrid(x_centers, y_centers)

    # surface
    ax.plot_surface(x_c,y_c,z_dem_2D, cmap = 'terrain', rstride=1, cstride=1, antialiased=True, shade=False)
            # x_centers,y_centers
    #ax.view_init(30, 60) # elev, azimuth, roll
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z') # [m agl.]')
    #ax.set_title('stacked '+args.variable)
    #ax.text2D(0.1,0.9,f"time={DumpInt*DumpGap}h",transform = ax.transAxes)

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.dem.3d.png', dpi = 300)

    # voxels
    #ax.voxels
    # https://stackoverflow.com/questions/73876939/plot-3d-grid-data-as-heat-map-using-matplotlib


def plot_3d_singlevardtsgl(vname,t):

    # cell centered variable
    data_3Dt = centered_var_generator(vname)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)
    
    pltsettings = variablescaling(data_3Dt[t], vname) # use 1 common scale of color for all layers and time steps
    varlabel = pltsettings[1]
    varrange = pltsettings[2]
    cmap = getattr(plt.cm, pltsettings[0])
    norm = pltsettings[3]

    print(data_3Dt.shape) # (93, 7, 100, 100)

    #x,y = np.meshgrid(np.arange(data_3Dt.shape[3]), np.arange(data_3Dt.shape[2]))
    x_c, y_c = np.meshgrid(x_centers, y_centers)

    for z in range(0,nz):
        
        data = data_3Dt[t,z]
        #print(data.shape)
        vmin = data.min()
        vmax = data.max()
        vmean = np.mean(data)
        print(f'layer {z} centered depth z=',z_centers[z],'m.agl; vmin,vmean,vmax',round(vmin,5),round(vmean,5),round(vmax,5))
    
        #if z in range(0,nz,1):
        ax.plot_surface(x_c,y_c,np.full_like(data, z_centers[z]), facecolors=cmap(data), norm=norm, rstride=1, cstride=1, antialiased=True, shade=False)
            # x_centers,y_centers
    #ax.view_init(30, 60) # elev, azimuth, roll
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z') # [m agl.]')
    #ax.set_title('stacked '+args.variable)
    #ax.text2D(0.1,0.9,f"time={DumpInt*DumpGap}h",transform = ax.transAxes)

    
    m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  #, norm=surf.norm)
    #vmin = data.min()
    #vmax = data.max()
    #ax.set_zlim(-1+DepthCenteredCellZ[0], 1+DepthCenteredCellZ[-1])
    m.set_clim(varrange[0],varrange[1])
    #m.set_clim(53.9982,53.9983)
    #m.set_clim(0,50)
    #print(vmin,vmax)
    #print(f'non all nul in cell : {np.where(data>-3.402823466385288e+38)}')
    plt.colorbar(m, ax=plt.gca(), label=varlabel)

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.{vname}.dt{t}.3d.png', dpi = 300)

def plot_proj2d_singlevardtslall(varname,projection, proj_idx_list, tmod, **kwargs):
    # plot for 3Dt variables only (ie. exclude overland_flow)
    vname = varname
    nbsl = len(proj_idx_list)
    
    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    nt = len(time_indexes)
    fig, axs = plt.subplots(nt,nbsl,figsize=(nbsl*6,nt*5),
                            gridspec_kw={'width_ratios': [1]*nbsl,
                                        'height_ratios': [1]*nt})

    projectioninfo = projscaling(vname, projection)
    x = projectioninfo[0]
    y = projectioninfo[1]
    xtit = projectioninfo[2]
    ytit = projectioninfo[3]
    #xlim = projectioninfo[4]
    #ylim = projectioninfo[5]

    if vname in pfb_out3dtvariables:
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    elif vname in pfb_cellcentered3dtvariables : 
        array_3Dt = centered_var_generator(vname)

    pltsettings = variablescaling(array_3Dt, vname) # use 1 common scale of color for all layers and time steps
    varcolor = pltsettings[0]
    varlabel = pltsettings[1]
    varrange = pltsettings[2]

    for idx in proj_idx_list:
        # print(idx, 'layer')
        array_2Dt = projected_array(array_3Dt,projection,idx)

        for t in range(nt):
            tidx = time_indexes[t]
            #print(t, 't axs loop', np.mean(array_2Dt[t]))
            im = axs[t,idx].pcolormesh(x, y, array_2Dt[tidx], shading='auto', cmap=varcolor, vmin = varrange[0], vmax=varrange[1] )

            if (projection == '2dxz') | (projection == '2dyz'):
                axs[t,idx].set_yticks(y)
            #axs[t,v].set_aspect(10000)
            elif projection == '2dxy':
                axs[t,idx].set_aspect('equal')
            
            axs[t,idx].set_xlabel(xtit)
            #axs[t,v].set_xlim(xlim[0],xlim[-1])
            axs[t,idx].set_ylabel(ytit)
            #axs[t,v].set_ylim(ylim[0],ylim[-1])

            realt = dt_real[tidx]
            axs[t,idx].set_title(f't={realt}h',loc='left') #, ha='center', va='center', transform=axs[t, v].transAxes)
            axs[t,idx].set_title(f'layer {idx}')
            fig.colorbar(im, ax=axs[t,idx], orientation='vertical') #, fraction=0.5, pad=0.04)
        
    #plt.text(f'{varlabel}[t,layer] for run {foldername}, in plan {projection}')
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.{vname}.dt{tmod}{nt}.{projection}.png', dpi = 300)
    plt.close()

def plot_proj2d_multivardtall(runname,projection, proj_idx, tmod, **kwargs):

    pfb_outvariables = pfb_out3dtvariables
    if (projection == '2dxy') & (proj_idx == nz-1):
        pfb_outvariables = pfb_out3dtvariables + pfb_out2dtvariables
    
    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    nt = len(time_indexes)
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
        #xlim = projectioninfo[4]
        #ylim = projectioninfo[5]

        if vname in pfb_out3dtvariables:
            if vname == 'vel_norm': 
                velx = readpfblist_to_3Dtarray(path, runname, 'velx')
                vely = readpfblist_to_3Dtarray(path, runname, 'vely')
                velz = readpfblist_to_3Dtarray(path, runname, 'velz')
                array_3Dt = vector_norm(velx,vely,velz)
            else:
                array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
            array_2Dt = projected_array(array_3Dt,projection,proj_idx)
        elif vname in pfb_out2dtvariables:
            array_2Dt = readpfblist_to_2Dtarray(path,runname,vname) #3Dt, future 2Dt
        #print(array_2Dt.shape)
        
        pltsettings = variablescaling(array_2Dt, vname)
        varcolor = pltsettings[0]
        varlabel = pltsettings[1]
        varrange = pltsettings[2]

        #array_fantom = np.empty(array_2Dt[0].shape)

        for t in range(nt):
            tidx = time_indexes[t]
            # print(t, f'{tidx} axs loop', np.mean(array_2Dt[tidx]))
            im = axs[t,v].pcolormesh(x, y, array_2Dt[tidx], shading='auto', cmap=varcolor, vmin = varrange[0], vmax=varrange[1] )

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
            axs[t, 0].set_title(f't={dt_real[tidx]}h',loc='left') #, ha='center', va='center', transform=axs[t, v].transAxes)
            fig.colorbar(im, ax=axs[t,v], orientation='vertical') #, fraction=0.5, pad=0.04)
            #fig.clim(varrange[0],varrange[1])

        axs[0,v].set_title(f'{varlabel}')
        
        #im_fantom = axs[t+1,v].pcolormesh(x,y,array_fantom)
        #axs[t+1,v].set_aspect('0.001')
        #axs[t+1,v].axis('off')
        #fig.colorbar(im, ax=axs[t+1,v], orientation='horizontal', label = varlabel, fraction=0.5, pad=0.04)
    
    #plt.title(f'Sequence of PF direct outputs for run {foldername}, in plan {projection}')
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.varall.dt{tmod}{nt}.{projection}{proj_idx}.png', dpi = 300)
    plt.close()

def plot_ZXoft(runname, Xcp, tmod, **kwargs):

    X = Xcp[0]
    nbX = X.shape[0]
    loc_cp = Xcp[1]
    position_layer = nz-1

    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    nt = len(time_indexes)
    tarray = dt_real[time_indexes]

    nbvar = len(pfb_out3dtvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    for v in range(nbvar):
        vname = pfb_out3dtvariables[v]
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        varlabel = pltsettings[1]

        for x in range(nbX):
            array_Xt = array_3Dt[time_indexes,position_layer,X[x,0],X[x,1]]
            axs[v].plot(tarray,array_Xt,label=loc_cp[x]+f'({position_layer},{X[x,0]},{X[x,1]})', marker='.')
    
        axs[v].grid(True)
        axs[v].legend()
        axs[v].set_xlabel('t [h]')
        axs[v].set_ylabel(f'{varlabel}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.varall.dt{tmod}.X.png', dpi = 300)

def plot_Zxoft(runname, Xcp, tmod,**kwargs):

    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    
    nbvar = len(pfb_out3dtvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    tarray = dt_real[time_indexes]

    for v in range(nbvar):
        vname = pfb_out3dtvariables[v]
        array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        varlabel = pltsettings[1]

        for z in range(nz):
            array_Xt = array_3Dt[time_indexes,z,X[0],X[1]]
            axs[v].plot(tarray,array_Xt,label=f'layer {z}', marker='.')
    
        axs[v].grid(True)
        axs[v].legend()
        axs[v].set_xlabel('t [h]')
        axs[v].set_ylabel(f'{varlabel}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.varall.dt{tmod}.{Xcp[1]}z.png', dpi = 300)

####################################
# CONVERGENCE PROXY

def plot_Hxoft(runname, Xcp, tmod,**kwargs):
    # compare pressure-assumed linear, vs, water table, (vs BC if point on y-lower)
    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    tarray = dt_real #np.linspace(0,nt-1,nt)
    
    nbvar = 1
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    nt = len(time_indexes)
    tarray = dt_real[time_indexes]
    #print(tarray)

    ## pseudo hydraulic head # h = p/(rho g) + z
    #vname = 'press'
    #array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    #pltsettings = variablescaling(array_3Dt, vname)
    #varlabel = pltsettings[1]

    #for z in range(nz):
    #    array_Xt = array_3Dt[:,z,X[0],X[1]] + z_centers[z] # h = p/(rho g) + z
    #    axs.plot(tarray,array_Xt,label=f'p(l={z})', marker='.')

    # hydraulic head # h = p/(rho g) + z + u^2/2
    array_3Dt = hydraulichead_3Dt()
    for z in range(nz):
        array_Xt = array_3Dt[time_indexes,z,X[0],X[1]]
        axs.plot(tarray,array_Xt,label=f'h(l={z})', marker='.')

    # Water table
    vname = 'water_table'
    v = 2
    array_3Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
    pltsettings = variablescaling(array_3Dt, vname)
    #varlabel = pltsettings[1]
    array_Xt = array_3Dt[time_indexes,X[0],X[1]]
    axs.plot(tarray,array_Xt,label='pf-wt', linestyle=':', color='k', marker='.')
    
    # PBC y-lower
    v = 3
    ylower = press_BC['y-lower__alltime'] * np.ones(nt)
    axs.plot(tarray,ylower,label='BCylow')

    axs.grid(True)
    axs.legend()
    axs.set_xlabel('t [h]')
    axs.set_ylabel('h [m agl]')

    plt.title('Hydraulic head') # and water table(<0)')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.Hxoft.dt{tmod}.{loc_cp}.png', dpi = 300)

def plot_HTofy(runname, tmod, **kwargs):
    # compare water table in space (hyp. steady state)    
    
    vname = 'water_table'
    array_2Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    #tarray = dt_real[time_indexes]
    
    # radius of influence - in steady state
    R = 1000
    nR = int(R/dx)

    fig, axs = plt.subplots(2,1)

    for t in time_indexes:
        array_1Dtfin = array_2Dt[t,:,int(nx/2)]
        axs[0].plot(y_centers[:nR],array_1Dtfin[:nR], label = f't={dt_real[t]}h', marker='.')
        axs[1].plot(y_centers,array_1Dtfin, label = f't={int(dt_real[t])}h')

    R_tot = dy*ny
    axs[0].plot(y_centers[:nR], unconfined_Dupuit(20/100,R_tot,1e-2,8,y_centers[:nR]), color ='k', linestyle =':', label='Dupuit eqt')
    axs[0].plot(y_centers[:nR], unconfined_PBC(8,7,R_tot,y_centers[:nR]), color ='b', linestyle =':', label='PBC eqt')
    axs[1].plot(y_centers, unconfined_Dupuit(20/100,R_tot,1e-2,8,y_centers), color ='k', linestyle =':', label='Dupuit eq')
    axs[1].plot(y_centers, unconfined_PBC(8,7,R_tot,y_centers), color ='b', linestyle =':', label='PBC eq')

    axs[0].grid(True)
    axs[1].grid(True)
    axs[0].set_aspect(100)
    axs[1].legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5)) #"Location","southeast")
    axs[0].set_xlabel('y [m]')
    axs[1].set_xlabel('y [m]')
    axs[0].set_ylabel('h [m agl]')
    axs[1].set_ylabel('h [m agl]')

    #plt.title(f'Water table')
    plt.savefig(f'{path_fig}{foldername}.HTofy.R{R_tot}.dt{tmod}.png', dpi = 300)


def plot_compared_surfaceflow(runname, Xcp, mod,**kwargs):
    # compare overland_flow (Kinematic wave equation) vs flow discharge in surface layer, at chosen point Xcp
    # ref to (Rousseau, 2020)

    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    tarray = dt_real #np.linspace(0,nt-1,nt)
    
    nbvar = 1
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))
    time_indexes = tidxscaling(dt_real,mod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    tarray = dt_real[time_indexes]

    # discharge
    # hypothesis that flow velocity in surface layer is the same as the one in the shear stress layer of free surface flow
    vname = 'discharge'
    array_3Dt = dischargenorm_3Dt()
    array_2Dt = projected_array(array_3Dt,"2dxy",nz-1)
    array_Xt = array_2Dt[time_indexes[:],X[0],X[1]]
    axs.plot(tarray,array_Xt, marker='.', label = r'$Q_{subsurface}$')

    # of
    vname = 'overland_flow'
    array_2Dt = readpfblist_to_2Dtarray(path,runname,vname)
    #pltsettings = variablescaling(array_2Dt, vname)
    #varlabel = pltsettings[1]

    array_Xt = array_2Dt[time_indexes,X[0],X[1]]
    axs.plot(tarray,array_Xt, marker='.',label=r'$Q_{overland flow}$')

    axs.grid(True)
    axs.legend()
    axs.set_xlabel('t [h]')
    axs.set_ylabel(r'$Q$ [m$^3$/h]')

    plt.title(f'Comparison of flux in surface layer and overland_flow in {loc_cp}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.fluxcomp.dt{mod}.{loc_cp}.png', dpi = 300)

def plot_boundarychecks(mod,**kwargs):
    # for all faces in press_BC list - for CV sum should go to 0 ('balanced model')
    # normal to face is directed outward of domain

    sum = np.zeros(nt)

    fig, ax = plt.subplots(1,1,figsize=(6,5))
    time_indexes = tidxscaling(dt_real,mod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    tarray = dt_real[time_indexes]

    # 'x-lower'
    array_Xt = (-1)*layer_sum1Dt(projected_array(discharge_3Dt('x'),'2dyz',0))
    sum = array_Xt
    #ax.plot(tarray,array_Xt[time_indexes],label='x-lower', marker='.')
    # 'x-upper'
    array_Xt = layer_sum1Dt(projected_array(discharge_3Dt('x'),'2dyz',nx))
    sum += array_Xt
    #ax.plot(tarray,array_Xt[time_indexes],label='x-upper', marker='.')
    # 'y-lower'
    array_Xt = (-1)*layer_sum1Dt(projected_array(discharge_3Dt('y'),'2dxz',0))
    sum += array_Xt
    ax.plot(tarray,array_Xt[time_indexes],label='y-lower', marker='.')
    # 'y-upper'
    array_Xt = layer_sum1Dt(projected_array(discharge_3Dt('y'),'2dxz',ny))
    sum += array_Xt
    ax.plot(tarray,array_Xt[time_indexes],label='y-upper', marker='.')
    # 'z-lower'
    array_Xt = (-1)*layer_sum1Dt(projected_array(discharge_3Dt('z'),'2dxy',0))
    sum += array_Xt
    #ax.plot(dt_real,array_Xt[time_indexes],label='z-lower', marker='.')
    # 'z-upper' - put on a secondary axis !
    array_Xt = layer_sum1Dt(projected_array(discharge_3Dt('z'),'2dxy',nz))
    sum += array_Xt
    ax.plot(tarray, array_Xt[time_indexes], label='z-upper', linestyle='', marker='+')
    
    # 'sum'
    ax.plot(tarray,sum[time_indexes],label='sum', color='k',linestyle=':')

    ax.grid(True)
    ax.legend()
    ax.set_xlabel('t [h]')
    ax.set_ylabel(r'$Q_{out}$ [m$^3$/h]')

    plt.title('Outflow of box-boundary faces')
    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.outflow.Bf-low.dt{mod}.png', dpi = 300)

def plot_MXFig3(): #runname, tmod, **kwargs): #MXFig3

    vname = 'water_table'
    array_2Dt = readpfblist_to_2Dtarray(path,runname,vname)
    #time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    time_indexes = [0,2,120]
    colors = ['k','b','r']
    i = 0

    fig, ax = plt.subplots(1,figsize=(5,5))

    for t in time_indexes:
        array_1Dtfin = array_2Dt[t,int(ny/2),:]
        ax.plot(x_centers,array_1Dtfin, label = f't={dt_real[t]}h', marker='.', color=colors[i])
        i +=1
    
    ax.grid(True)
    ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    ax.set_xlabel('x [m]')
    ax.set_ylabel('wt [m]')

    #plt.title(f'Water table')
    plt.savefig(f'{path_fig}{foldername}.MXFig3.png', dpi = 300)


def plot_HTofx(): #runname, tmod, **kwargs): #MXFig3

    time_indexes = [0,2,120]

    vname = 'water_table'
    array_2Dt = readpfblist_to_2Dtarray(path,runname,vname)
    #time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    colors = ['k','b','r']
    c = 0

    array_3Dt = hydraulichead_3Dt()
    rgz = [0, int(nz/2), nz-1]

    fig, ax = plt.subplots(1,figsize=(5,5))

    for t in time_indexes:
        array_WT1Dt = -1 * array_2Dt[t,int(ny/2),:]
        ax.plot(x_centers,array_WT1Dt, color = colors[c],label = f't={dt_real[t]}h')
        array_H1Dt = array_3Dt[t,rgz,int(ny/2),:]
        ax.plot(x_centers, array_H1Dt[0], linestyle = 'dotted', color = colors[c])
        ax.plot(x_centers, array_H1Dt[1], linestyle = 'dashed', color = colors[c])
        ax.plot(x_centers, array_H1Dt[2], linestyle = 'dashdot', color = colors[c])
        c +=1
    
    ax.grid(True)
    ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    ax.set_xlabel('x [m]')
    ax.set_ylabel('h [m as]')

    #plt.title(f'Water table')
    plt.savefig(f'{path_fig}{foldername}.HTofx.dtFig3.png', dpi = 300)


def plot_HTofybis(): #runname, tmod, **kwargs): #MXFig3

    time_indexes = [0,4,nt-1]

    vname = 'water_table'
    array_2Dt = readpfblist_to_2Dtarray(path,runname,vname)
    #time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    colors = ['k','b','r']
    c = 0

    array_3Dt = hydraulichead_3Dt()
    rgz = [0, int(nz/2), nz-1]

    fig, ax = plt.subplots(1,figsize=(5,5))

    for t in time_indexes:
        array_WT1Dt = -1 * array_2Dt[t,:,int(nx/2)]
        ax.plot(y_centers,array_WT1Dt, color = colors[c],label = f't={dt_real[t]}h')
        array_H1Dt = array_3Dt[t,rgz,:,int(nx/2)]
        ax.plot(y_centers, array_H1Dt[0], linestyle = 'dotted', color = colors[c])
        ax.plot(y_centers, array_H1Dt[1], linestyle = 'dashed', color = colors[c])
        ax.plot(y_centers, array_H1Dt[2], linestyle = 'dashdot', color = colors[c])
        c +=1
    
    ax.plot(y_centers, unconfined_PBC(8,7,ny*dy,y_centers), color ='g', label='ss.eq')

    ax.grid(True)
    ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    ax.set_xlabel('x [m]')
    ax.set_ylabel('h [m as]')

    #plt.title(f'Water table')
    plt.savefig(f'{path_fig}{foldername}.HTofy.dtFig3.png', dpi = 300)

################################################################
# THEORETICAL EQUATION

def unconfined_Dupuit(Q,R,K,ho,y_array): # steady stata
    return (Q/(K*np.pi)*np.log10(y_array/R) + ho**2)**(1/2) + z_faces[0]

def unconfined_PBC(hup,hlow,Ly,y): #steady state
    return ((hup**2-hlow**2)/Ly*y+hlow**2)**(1/2) + z_faces[0]

################################################################

# Nomenclature of variables

pfb_out3dtvariables = ['press','satur','velx','vely','velz']
pfb_cellcentered3dtvariables = ['press', 'satur', 'vx_c','vy_c','vz_c','v', 'H'] #'subsurface_storage' - H for hydraulic head
pfb_out2dtvariables = ['water_table', 'overland_flow'] #'surface_storage', 'surface_storage', 'overland_bc_flux', 'overlandsum',
faces = ['x-lower', 'x-upper', 'y-lower', 'y-upper', 'z-lower', 'z-upper']
plot_projections = ['2dxy','2dxz','2dyz','1d','3d']

zidx = nz-1
yidx = 0
xidx = int(nx/2)

plot_3d_geom()

"""
#tmod = "all"
#tmod = "beg" # beginning
tmod = "int" # interval
tmod = "ed"
ndt = 3
tinit = dt_real[1] # initial time in hours - for plot
tfin = 120 #dt_real[-1]-1 # final time in hours - for plot
plot_proj2d_multivardtall(runname,'2dxz',yidx, tmod, c=ndt, d=tinit, e=tfin)
## plot_proj2d_multivardtall(runname,'2dxz',yidx, "all")

plot_proj2d_multivardtall(runname,'2dxy',zidx, tmod, c=ndt, d=tinit, e=tfin)
plot_proj2d_multivardtall(runname,'2dyz',xidx, tmod, c=ndt, d=tinit, e=tfin)

#plot_3d_singlevardtsgl('H', 7) # for cellcentered3dtvariables

#layerstoplot = [10,9,8,7,6,5,4,3,2,1,0]  #VM
layerstoplot = [6,5,4,3,2,1,0]
#plot_proj2d_singlevardtslall('H','2dxy',layerstoplot, tmod, c=ndt, d=tinit, e=tfin)


## VM check points
f_cp = '/home/patras/PF-Valmalenco/data/controlpoints.txt'
Xidx_cp = read_cpcsv(f_cp)
#print(Xidx_cp)
XP_ylower = [Xidx_cp[0][1],[Xidx_cp[1][1]]]
XP_center = [np.array([int(ny/2),int(nx/2)]),f'[{int(ny/2)},{int(nx/2)}]']

## 'Universal' check point
Xidx_cp = [np.array([[0,int(nx/2)],[1,int(nx/2)],[2,int(nx/2)],[int(ny/2),int(nx/2)],[ny-1,int(nx/5)],[ny-1,int(nx/2)]]),['P1','P2','P3','P4','P5','P6']]
XP_center = [np.array([int(ny/2),int(nx/2)]),f'[{int(ny/2)},{int(nx/2)}]']
XP_ylower = [np.array([0,int(nx/2)]),f'[0,{int(nx/2)}]']
XP_ylowerbis = [np.array([1,int(nx/2)]),f'[1,{int(nx/2)}]']

# Find max index in overland_flow # [[3, 67], 'X_ofmax'] for VM_V52
# Q_max(x=16750) = 1176 m^3/h = 0.32 m^3/s
array_forXPmax = readpfblist_to_2Dtarray(path,runname,'overland_flow')
max_of = array_forXPmax.max()
#print(max_of)
idx_max = np.where(array_forXPmax == array_forXPmax.max())
idx_np = np.array(idx_max)
idxlist = [arr[0] for arr in idx_np] # Convert the numpy arrays into a list of values
XP_maxof = [idxlist[1:],'X_ofmax']
#print(XP_maxof)


tmod = "all"
#tmod = "beg" # beginning
#tmod = "intall" # interval
#tinit = 0 # initial time in hours - for plot
#tfin = 5 # final time in hours - for plot
plot_ZXoft(runname,Xidx_cp,tmod, d=tinit, e=tfin) # X : multiple points
plot_Zxoft(runname,XP_ylower,tmod, d=tinit, e=tfin) # x : single point
plot_Zxoft(runname,XP_center,tmod, d=tinit, e=tfin)
plot_Hxoft(runname,XP_ylower,tmod, d=tinit, e=tfin) # Hxoft = H(x)(t), evolution of H(x) in time
plot_Hxoft(runname,XP_ylowerbis,tmod, d=tinit, e=tfin)
plot_Hxoft(runname,XP_center,tmod, d=tinit, e=tfin)
plot_compared_surfaceflow(runname,XP_maxof,tmod, d=tinit, e=tfin)
plot_boundarychecks(tmod,c=ndt, d=tinit, e=tfin)


#tmod = "fin"
ndt = 5
tmod = "int" # interval
tinit = dt_real[1] # initial time in hours - for plot
tfin = 120 #dt_real[-1]-1 # 5e7 # final time in hours - for plot #print(dt_real[-1])
plot_HTofy(runname,tmod, c=ndt, d=tinit, e=tfin ) # H(t) distributed in space along y axis (x=nx/2)
"""

#plot_MXFig3()
plot_HTofybis()
