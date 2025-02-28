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

from crs_to_idx import *
from read_log import *


# coment if not defined
plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
colorspalette = ['tab:grey','tab:blue','tab:cyan','tab:green','tab:olive','tab:orange','tab:brown','tab:red','tab:pink','tab:purple']

#####################################################################################"
# OUTPUT PATH

#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant

###############################################################################################
# INPUT PATH
# define the complete path to your output .pfidb and .pfb files
# such that <runname>.pfidb can be found in : f'{path}{foldername}{runname}.pfidb'

## DS
#path = '/home/patras/PF-Test/DumbSquare/outputs/'
#foldername = 'DSc1000z10s0.MX.DP23.IN0dt01_v55' #'DS.c100s1_v28'
#runname = 'DSc1000z10s0'

## MX
#path = '/home/patras/PF-Test/Maxwell2013/outputs/'
#foldername = 'MX.c1s1y3_v6' #'MX.c100s1bcx_v5'
#runname = 'MX.c1s1'

## VM
path = '/home/patras/PF-Valmalenco/outputs/BowlBox_UXZ/'
foldername = 'BowlBox_UXZ_TgwDpRFe-4K36e-10L1'
runname = 'PLT_box'

#path = '/home/patras/PF-Valmalenco/outputs/BowlBox_UXZ/'
#foldername = 'BowlBox_UXZ_K36e-3L6d' #'BowlBox_UXZ_TgwDpRFe-4K36e-10L1'
#runname = 'BowlBox_UXZ' #'BowlBox_UXZ' #'PLT_box'

#path = '/home/patras/PF-Valmalenco/outputs/'
#foldername = 'SeepBox_KSsgL11IC-0125' #'BowlBox_UXZ_K36e-3L6d' #'BowlBox_UXZ_TgwDpRFe-4K36e-10L1'
#runname = 'SeepBox' #'BowlBox_UXZ' #'PLT_box'

#path = '/home/patras/PF-Valmalenco/outputs/'
#foldername = 'BB.c500.VMRF0IC-0375' #'BowlBox_UXZ_K36e-3L6d' #'BowlBox_UXZ_TgwDpRFe-4K36e-10L1'
#runname = 'BB.c500.VMRF0IC-0375' #'BowlBox_UXZ' #'PLT_box'

#path = '/home/patras/PF-Valmalenco/outputs/'
#foldername = 'CLM_V52'
#runname = 'CLM_V5'

## VM
#path = '/home/patras/PF-Valmalenco/outputs/'
#foldername = 'RainRec_RFe-4C2-2IC-35'
#runname = 'RainRec_RFe-4C2-2IC-35'

## parflow test cases
#path = '/home/patras/PF-Test/LW/outputs/'
#foldername = 'LW_var_dz_spinup'
#runname = 'LW_var_dz_spinup'

# NB: Possible error message on Boundary Condition Patches, if so comment line 183,184.

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

X_origin = [280,0,0] # change with real CRS coordinate
print('Domain origin, south-west:', X_origin[0])

# domain rect cells
x_centers = dx*(np.linspace(1,nx,nx)-1/2)
y_centers = dy*(np.linspace(1,ny,ny)-1/2)
x_faces = dx*(np.linspace(0,nx,nx+1))
y_faces = dy*(np.linspace(0,ny,ny+1))
#print(x_centers.shape)
#print(y_centers.shape)

dt_real = dump_to_simulatedtimes_equivalent(path,foldername,runname)
print(f'dumptimes real {dt_real}')
Deltat = dump_timesteps(path,foldername,runname)
nt = len(dt_real)

#def layers_centereddepth():
dzscaled = np.concatenate((dz,[0]))
#print(dz)
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
print('abs slope_x : ', np.min(abs(slopex)), np.mean(abs(slopex)), np.max(abs(slopex)))
print('abs slope_y : ', np.min(abs(slopey)), np.mean(abs(slopey)), np.max(abs(slopey)))

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
print('z min, max', np.min(z_centers_asl_3D), np.max(z_centers_asl_3D))

dt_real = dump_to_simulatedtimes_equivalent(path,foldername,runname)
#print(f'dumptimes {dt_real}')
nt = len(dt_real)

#################################
# DISTRIBUTED INPUT PARAMETERS

pfbinput_variables = ['alpha','n','porosity','specific_storage','mannings','permx','permy','permz']
pfbinput_MRC = ['alpha','n','sres','ssat']

def readpfbinput_to_array(path,runname,variable):
    # elaborate 2D or 3D array for input values  pfbinput_variables
    if variable in pfbinput_MRC:
        filename = path + foldername + '/' + runname + ".out." + variable +".pfb"
        variable_array = read_pfb(filename)
    else:
        print("Error: not an input variable")
    return variable_array

porosity = data.computed_porosity
specific_storage = data.specific_storage
mannings = data.mannings #2D
permx = data.computed_permeability_x
permy = data.computed_permeability_y
permz = data.computed_permeability_z
perm_avg = np.mean(permz)
print(perm_avg)

alpha = readpfbinput_to_array(path,runname,'alpha')
n = readpfbinput_to_array(path,runname,'n')
sres = readpfbinput_to_array(path,runname,'sres')
ssat = readpfbinput_to_array(path,runname,'ssat')

#press_BC = data.pressure_boundary_conditions
#print('Pressure head boundary conditions : ', press_BC)
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
        variable2Dt = variable2Dt #*1/3600 # m^3/s

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
  
def bernoullihead_3Dt():
    # centered variable
    # nb. Bernoulli assumptions is not verified (not Newtonian), thus h is an estimate of the hydraulic head
    vel = velocitynorm_3Dt()
    press = readpfblist_to_3Dtarray(path, runname, 'press') # pressure head
    # H = p/rhog + z + u^2/rho2
    # Element-wise operation
    H = press + vel**2 + z_centers_3D
    return H

def hydraulichead_3Dt():
    # centered variable
    # nb. Bernoulli assumptions is not verified (not Newtonian), thus h is an estimate of the hydraulic head
    vel = velocitynorm_3Dt()
    press = readpfblist_to_3Dtarray(path, runname, 'press') # pressure head
    # h = p/rhog + z 
    # Element-wise operation
    h = press + vel**2 # + z_centers_3D
    return h

def hydraulichhead_3Dt_bis(): # BETTER USE hydraulichead_3Dt for H, this compute h 
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
        array3Dt = bernoullihead_3Dt()
    elif vname == 'h':
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

def diff_array2D(array1,array2):
    diff = array2 - array1
    meandiff = np.mean(diff)
    mindiff = np.min(diff)
    #minidx = np.where(diff==mindiff)
    maxdiff = np.max(diff)
    #minatt = [mindiff,minidx]
    #maxidx = np.where(diff==maxdiff)
    #maxatt = [maxdiff,maxidx]
    return diff,meandiff,mindiff,maxdiff

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
    arrayshape = array_2Dt.shape
    array_1Dt = np.zeros(arrayshape[0])
    for t in range(arrayshape[0]):
        array_1Dt[t] = np.mean(array_2Dt[t])
    return array_1Dt

def layer_mean3Dtto1Dt(array_3Dt):

    s = array_3Dt.shape()
    distt = s[0]

def layer_sum1Dt(array_2Dt):
    # average a 2D array for each timestep
    array_1Dt = np.zeros(nt)
    for t in range(nt):
        array_1Dt[t] = np.sum(array_2Dt[t])
    return array_1Dt

def maxinspace_3Dt(array_3Dt):

    array_1Dt = np.max(array_3Dt,3)
    array_1Dt = np.max(array_1Dt,2)
    array_1Dt = np.max(array_1Dt,1)

    return array_1Dt

def derivative3D_dvarofdt(array_3Dt):

    derivative_3Dt = np.empty((nt-1,nz,ny,nx))
    for t in range(1,nt):
        derivative_3Dt[t-1] = (array_3Dt[t] - array_3Dt[t-1])/Deltat[t] # forward derivative
    return derivative_3Dt

def derivative2D_dvarofdt(array_2Dt):

    derivative_2Dt = np.empty((nt-1,ny,nx))
    for t in range(1,nt):
        derivative_2Dt[t-1] = (array_2Dt[t] - array_2Dt[t-1])/Deltat[t] # forward derivative
    return derivative_2Dt

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
        varrange = [max(array.min(),-10), min(array.max(),10)] #min(array.max(),2)]
        #varrange = [-1.5,1.5] # [0.7, 0.1]
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
        varrange = [0, 600] #[array.min(), array.max()]
    elif variable == 'H':
        colorofmap = 'Blues'
        varlabel = 'H [m asl.]'
        varrange = [array.min(),array.max()]
        #varrange = [-1.0,1.0]
    elif variable == 'h':
        colorofmap = 'Blues'
        varlabel = 'h [m agl.]'
        if array.min()<0:
            normed_range = min(abs(array.min()), array.max())
            varrange = [-1*normed_range, normed_range]
        else:
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

#### FOR INPUTS (3 first function) AND OUTPUTS

def plot_3d_geom():

    pltsettings = variablescaling(z_dem_2D, 'DEM') # use 1 common scale of color for all layers and time steps
    varcolor = pltsettings[0]
    varlabel = pltsettings[1]
    varrange = pltsettings[2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)
    
    #x,y = np.meshgrid(np.arange(data_3Dt.shape[3]), np.arange(data_3Dt.shape[2]))
    x_c, y_c = np.meshgrid(x_centers, y_centers)

    # surface
    ax.plot_surface(x_c,y_c,z_dem_2D, cmap = varcolor, rstride=1, cstride=1, antialiased=True, shade=False, vmin = varrange[0], vmax=varrange[1])
            # x_centers,y_centers
    ax.view_init(15, -100) # elev, azimuth, roll
    
    ax.set_box_aspect(aspect=(2, 2, 1), zoom=1.3)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z [m asl.]')
    #ax.set_title('stacked '+args.variable)
    #ax.text2D(0.1,0.9,f"time={DumpInt*DumpGap}h",transform = ax.transAxes)

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.dem.3d.png', dpi = 300)

    # voxels
    #ax.voxels
    # https://stackoverflow.com/questions/73876939/plot-3d-grid-data-as-heat-map-using-matplotlib

plot_3d_geom()

def plot_proj2d_singlevar_geom():
    # plot for 3Dt variables only (ie. exclude overland_flow)
    
    # dem, slopex, slopey
    fig, axs = plt.subplots(1,3,figsize=(3*6,5)) #

    projectioninfo = projscaling('press', '2dxy')
    x = projectioninfo[0]
    y = projectioninfo[1]
    xtit = projectioninfo[2]
    ytit = projectioninfo[3]
    #xlim = projectioninfo[4]
    #ylim = projectioninfo[5]

    # DEM
    pltsettings = variablescaling(z_dem_2D, 'DEM') # use 1 common scale of color for all layers and time steps
    varcolor = pltsettings[0]
    varlabel = pltsettings[1]
    varrange = pltsettings[2]

    im = axs[0].pcolormesh(x, y, z_dem_2D, shading='auto', cmap= varcolor, vmin = varrange[0], vmax=varrange[1] )
    #axs[0,1].set_aspect('equal')      
    axs[0].set_xlabel(xtit)
    #axs[t,v].set_xlim(xlim[0],xlim[-1])
    axs[0].set_ylabel(ytit)
    #axs[t,v].set_ylim(ylim[0],ylim[-1])
    axs[0].set_title(varlabel)
    fig.colorbar(im, ax=axs[0], orientation='vertical') #, fraction=0.5, pad=0.04)

    # SLOPE X
    pltsettings = variablescaling(slopex[0], 'slopex') # use 1 common scale of color for all layers and time steps
    varcolor = pltsettings[0]
    varlabel = pltsettings[1]
    varrange = pltsettings[2]

    im = axs[1].pcolormesh(x, y, slopex[0], shading='auto', cmap=varcolor, vmin = varrange[0], vmax=varrange[1] )
    #axs[1].set_aspect('equal')      
    axs[1].set_xlabel(xtit)
    #axs[t,v].set_xlim(xlim[0],xlim[-1])
    axs[1].set_ylabel(ytit)
    #axs[t,v].set_ylim(ylim[0],ylim[-1])
    axs[1].set_title(varlabel)
    fig.colorbar(im, ax=axs[1], orientation='vertical') #, fraction=0.5, pad=0.04)

    # SLOPE Y
    pltsettings = variablescaling(slopey[0], 'slopey') # use 1 common scale of color for all layers and time steps
    varcolor = pltsettings[0]
    varlabel = pltsettings[1]
    varrange = pltsettings[2]

    im = axs[2].pcolormesh(x, y, slopey[0], shading='auto', cmap=varcolor, vmin = varrange[0], vmax=varrange[1] )
    #axs[0,2].set_aspect('equal')      
    axs[2].set_xlabel(xtit)
    #axs[t,v].set_xlim(xlim[0],xlim[-1])
    axs[2].set_ylabel(ytit)
    #axs[t,v].set_ylim(ylim[0],ylim[-1])
    axs[2].set_title(varlabel)
    fig.colorbar(im, ax=axs[2], orientation='vertical') #, fraction=0.5, pad=0.04)

    #plt.text(f'{varlabel}[t,layer] for run {foldername}, in plan {projection}')
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.topo.2dxy.png', dpi = 300)
    plt.close()

plot_proj2d_singlevar_geom()

def plot_proj2d_MRC():
    # plot for 3Dt variables only (ie. exclude overland_flow)
    
    sl = nz
    MRC_var = len(pfbinput_MRC)
    fig, axs = plt.subplots(MRC_var,sl,figsize=(sl*6,MRC_var*5))

    # cell-centered variables
    projectioninfo = projscaling('press', '2dxy')
    x = projectioninfo[0]
    y = projectioninfo[1]
    xtit = projectioninfo[2]
    ytit = projectioninfo[3]
    #xlim = projectioninfo[4]
    #ylim = projectioninfo[5]

    for i in range(MRC_var):
        vname = pfbinput_MRC[i]
        varray = readpfbinput_to_array(path,runname,vname)
        pltsettings = variablescaling(varray, vname)
        varcolor = pltsettings[0]
        varlabel = pltsettings[1]
        varrange = pltsettings[2]

        for j in range(sl):
            im = axs[i,j].pcolormesh(x, y, varray[j], shading='auto', cmap= varcolor, vmin = varrange[0], vmax=varrange[1] )
            #axs[0].set_xlabel(xtit)
            #axs[0].set_ylabel(ytit)
            #axs[t,v].set_ylim(ylim[0],ylim[-1])
            axs[i,j].set_title(varlabel)
        fig.colorbar(im, ax=axs[i,j], orientation='vertical') #, fraction=0.5, pad=0.04)

    #plt.text(f'{varlabel}[t,layer] for run {foldername}, in plan {projection}')
    plt.tight_layout()
    plt.savefig(f'{path_fig}{foldername}.mrc.2dxyzall.png', dpi = 300)
    plt.close()

#plot_proj2d_MRC()

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
        pfb_outvariables = pfb_out3dtvariables[:2] + [pfb_out2dtvariables[0]] + pfb_out3dtvariables[2:] + [pfb_out2dtvariables[1]]
    
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
        
        pltsettings = variablescaling(array_2Dt[time_indexes], vname)
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
            axs[t, v].set_title(f't={dt_real[tidx]}h',loc='left') #, ha='center', va='center', transform=axs[t, v].transAxes)
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
    plt.close()

def plot_Zxoft(runname, Xcp, tmod,**kwargs):

    X = Xcp[0] # single point
    loc_cp = Xcp[1]
    
    nbvar = len(pfb_out3dtvariables)
    fig, axs = plt.subplots(1,nbvar,figsize=(nbvar*6,5))

    tfin = kwargs.get('e',float)
    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=tfin)
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
    plt.savefig(f'{path_fig}{foldername}.varall.dt{tmod}{tfin}.{Xcp[1]}z.png', dpi = 300)
    plt.close()

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
    array_3Dt = bernoullihead_3Dt()
    pltsettings = variablescaling(array_3Dt, 'H')
    varlabel = pltsettings[1]

    for z in range(nz):
        array_Xt = array_3Dt[time_indexes,z,X[0],X[1]]
        axs.plot(tarray,array_Xt,label=f'layer {z}', marker='.')
        
    if np.mean(array_Xt)<0:
        # Water table
        vname = 'water_table'
        v = 2
        array_3Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
        pltsettings = variablescaling(array_3Dt, vname)
        #varlabel = pltsettings[1]
        array_Xt = array_3Dt[time_indexes,X[0],X[1]]
        axs.plot(tarray,array_Xt,label='pf-wt', linestyle=':', color='k', marker='.')
    
    # PBC y-lower
    #v = 3
    #ylower = press_BC['y-lower__alltime'] * np.ones(nt)
    #axs.plot(tarray,ylower,label='BCylow')

    axs.grid(True)
    axs.legend()
    axs.set_xlabel('t [h]')
    axs.set_ylabel(varlabel)

    #plt.title('Hydraulic head') # and water table(<0)')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.Hxoft.dt{tmod}.{loc_cp}.png', dpi = 300)
    plt.close()

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
    plt.close()

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
    array_2Dt_l10 = projected_array(array_3Dt,"2dxy",nz-1)
    array_2Dt_l9 = projected_array(array_3Dt,"2dxy",nz-2)
    array_2Dt_l8 = projected_array(array_3Dt,"2dxy",nz-3)
    array_2Dt_l7 = projected_array(array_3Dt,"2dxy",nz-4)
    array_2Dt_l6 = projected_array(array_3Dt,"2dxy",nz-5)

    array_Xt_l10 = array_2Dt_l10[time_indexes[:],X[0],X[1]]
    array_Xt_l9 = array_2Dt_l9[time_indexes[:],X[0],X[1]]
    array_Xt_l8 = array_2Dt_l8[time_indexes[:],X[0],X[1]]
    array_Xt_l7 = array_2Dt_l7[time_indexes[:],X[0],X[1]]
    array_Xt_l6 = array_2Dt_l6[time_indexes[:],X[0],X[1]]
    
    axs.plot(tarray,array_Xt_l10, marker='.', label = r'$Q_{subsurface-l10}$')
    axs.plot(tarray,array_Xt_l9, marker='.', label = r'$Q_{subsurface-l9}$')
    axs.plot(tarray,array_Xt_l8, marker='.', label = r'$Q_{subsurface-l8}$')
    axs.plot(tarray,array_Xt_l7, marker='.', label = r'$Q_{subsurface-l7}$')
    axs.plot(tarray,array_Xt_l6, marker='.', label = r'$Q_{subsurface-l6}$')

    array_Qsumt_l106 = array_Xt_l10 + array_Xt_l9 + array_Xt_l8 + array_Xt_l7 + array_Xt_l6
    axs.plot(tarray,array_Qsumt_l106, marker='.', label = r'$Q_{subsurface-\Sigma l10-6}$')

    # of
    vname = 'overland_flow'
    array_2Dt = readpfblist_to_2Dtarray(path,runname,vname)
    #pltsettings = variablescaling(array_2Dt, vname)
    #varlabel = pltsettings[1]

    array_Xt = array_2Dt[time_indexes,X[0],X[1]]
    axs.plot(tarray,array_Xt, marker='.', color='k', linestyle=':', label=r'$Q_{overland flow}$')

    axs.grid(True)
    axs.legend(loc='upper right')
    axs.set_xlabel('t [h]')
    axs.set_ylabel(r'$Q$ [m$^3$/h]')
    axs.set_ylim(0,250)

    #plt.title(f'Comparison of flux in surface layer and overland_flow in {loc_cp}')

    #plt.show()
    plt.savefig(f'{path_fig}{foldername}.fluxcomp.dt{mod}.{loc_cp}.png', dpi = 300)
    plt.close()

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
    plt.close()

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

    array_3Dt = bernoullihead_3Dt()
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
    plt.close()

def plot_HTofybis(): #runname, tmod, **kwargs): #MXFig3

    time_indexes = [0,100,nt-1] #nt-1

    vname = 'water_table'
    array_2Dt = readpfblist_to_2Dtarray(path,runname,vname)
    #time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    colors = ['k','b','r']
    c = 0

    array_3Dt = bernoullihead_3Dt()
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


def cvproxy_toknownsolution(tmod,cvratio,**kwargs):
    # cv proxy that evolves toward 0 when simu CV to analytical/known solution
    
    vname = 'water_table'
    array_2Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
    tmod = 'all'
    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    #tarray = dt_real[time_indexes]
    array_known = unconfined_PBC_2D()

    diff_2Dt = np.ones((nt,ny,nx))
    cvproxy = np.ones((nt,3))
    for t in time_indexes:
        diff_res = diff_array2D(array_2Dt[t],array_known)
        diff_2Dt[t] = diff_res[0]/array_known
        cvproxy[t]=[np.mean(diff_2Dt[t]),np.mean(diff_2Dt[t,:,int(nx/2)]),np.max(diff_2Dt[t,:,int(nx/2)])]

    tcv = np.where(cvproxy[:,1]<cvratio)
    if len(tcv[0])>0:
        print(f'time of cv (<{cvratio}) t[{tcv[0][0]}] = {dt_real[tcv[0][0]]}')

    fig, axs = plt.subplots(1)
    axs.plot(dt_real,cvproxy[:,0])
    axs.plot(dt_real,cvproxy[:,1])
    axs.plot(dt_real,cvproxy[:,2])
    axs.plot(y_centers,cvratio*np.ones(y_centers.shape),color='k',linestyle=':',label=r'5$\%$')

    axs.grid(True)
    #axs.legend(loc='lower right') #best', bbox_to_anchor=(0.5, 0., 0.5, 0.5)) #"Location","southeast")
    axs.set_xlabel('t [h]')
    axs.set_ylabel(r'$\bar{\Delta h/h}$ [m agl]')

    #plt.title(f'Water table')
    plt.savefig(f'{path_fig}{foldername}.DHToft.m.png', dpi = 300)
    plt.close()

    return tcv

def plot_DHTofy_toknownsolution(tmod,**kwargs):
    # cv proxy that evolves toward 0 when simu CV to analytical/known solution
    
    vname = 'water_table'
    array_2Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    #tarray = dt_real[time_indexes]
    
    array_known = unconfined_PBC_2D()
    # radius of influence - in steady state

    fig, axs = plt.subplots(3,1,figsize=(1*5,3*2))

    for t in time_indexes:
        print(t, 'time idx')
        array_2D = array_2Dt[t,:,:]
        diffres = diff_array2D(array_2D,array_known)
        diff_2D = diffres[0]
        diff_1D = diff_2D[0,:,int(nx/2)]
        diffratio_1D = diff_1D/array_known[0,:,int(nx/2)]
        print(diffres[1:])
        axs[1].plot(y_centers,diff_1D)
        axs[0].plot(y_centers,array_2D[:,int(nx/2)], label = f't={round(dt_real[t],0)}h') #, marker='.')
        axs[2].plot(y_centers,diffratio_1D)
    
    axs[2].plot(y_centers,cv_ratio*np.ones(y_centers.shape),color='k',linestyle=':',label=fr'{cv_ratio}$\%$')

    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)
    axs[0].legend(loc='lower right') #best', bbox_to_anchor=(0.5, 0., 0.5, 0.5)) #"Location","southeast")
    axs[2].legend(loc='lower right') #best', bbox_to_anchor=(0.5, 0., 0.5, 0.5)) #"Location","southeast")
    axs[2].set_xlabel('y [m]')
    axs[1].set_ylabel(r'$\Delta$h [m agl]')
    axs[0].set_ylabel(r'h [m agl]')
    axs[2].set_ylabel(r'$\Delta$h/h [m agl]')

    #plt.title(f'Water table')
    plt.savefig(f'{path_fig}{foldername}.DHTofy.dt{tmod}.png', dpi = 300)
    plt.close()

def plot_wtderivative():

    print('Derivative evaluation')
    vname = 'water_table'
    array_2Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)
    # temporal derivative of wt
    deriv_wt = derivative2D_dvarofdt(array_2Dt)
    # mean in space of temporal derivative of wt
    mean_deriv_wt = layer_mean1Dt(deriv_wt)
    print('initial derivative', mean_deriv_wt[0])
    print('final derivative',mean_deriv_wt[-1])
    abs_mean_deriv_wt = abs(mean_deriv_wt)

    threshold = 1e-8
    idxt = np.where(abs_mean_deriv_wt<threshold)
    if len(idxt[0]>0):
        idx = idxt[0][0]
        print('idx lower 1e-8', idx)
    else:
        idx = nt-2
        print('threshold 1e-8 not reached')
    
    if runname == 'DSc1000z10s0':
        cv_ratio = 0.1
        tcv = cvproxy_toknownsolution("all",cv_ratio)
        if len(tcv[0]>0):
        #    print(f'time of cv (<5%) t[{tcv[0][0]}] = {dt_real[tcv[0][0]]}')
            idx_ks = tcv[0][0]
        #else:
        #    idx = nt-2

    fig, ax = plt.subplots(1)
    ax.plot(dt_real[1:idx+2], abs_mean_deriv_wt[:idx+1], label = r'dh/dt$_{mean,X}$', marker='.')
    # ax.plot(dt_real[1:idx+2], 1e-8*np.ones((idx+1)), label = r'threshold 1e-6', marker='.')
    if (runname == 'DSc1000z10s0') and len(tcv[0]>0):
        #ax.plot(dt_real[1:idx+2], abs_mean_deriv_wt[:idx+1], label = r'dh/dt$_{mean,X}$', marker='.')
        ax.plot(dt_real[1:idx+2], abs_mean_deriv_wt[idx_ks]*np.ones((idx+1)), label = f'dh/dt(t|dh/h<{cv_ratio})') #, marker='.')
    #else:
    #    ax.plot(dt_real[1:idx+2], abs_mean_deriv_wt[:idx+2], label = r'dh/dt$_{mean,X}$', marker='.')

    ax.grid(True)
    ax.legend(loc='upper right') #, bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    plt.yscale("log")
    ax.tick_params(axis='y', labelrotation=90)
    ax.set_xlabel('t [h]')
    ax.set_ylabel(r'dh/dt$_{mean,X}$')

    plt.savefig(f'{path_fig}{foldername}.dhdt-mean.png', dpi = 300)
    plt.close()

def plot_YofX(Yvarname, Xvarname, Xcp): # hysteretic curve

    #Yvarname = 'H'
    #Xvarname = 'satur'

    X = Xcp[0] # single point
    loc_cp = Xcp[1]

    vname = Yvarname
    if vname in pfb_cellcentered3dtvariables:
        Y_array_3Dt = centered_var_generator(vname)
    else:
        Y_array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    Y_array_Xt = Y_array_3Dt[:,:,X[0],X[1]]
    pltsettings = variablescaling(Y_array_3Dt, vname)
    Y_varlabel = pltsettings[1]
    # X
    vname = Xvarname
    X_array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    X_array_Xt = X_array_3Dt[:,:,X[0],X[1]]
    pltsettings = variablescaling(X_array_3Dt, vname)
    X_varlabel = pltsettings[1]

    fig, ax = plt.subplots()
    for z in range(nz-1):
        ax.plot(X_array_Xt[:,z],Y_array_Xt[:,z],label=f'layer {z}', marker='.', color=colorspalette[z])
        ax.plot(X_array_Xt[0,z],Y_array_Xt[0,z], marker='+', markersize=6, color=colorspalette[z])

    ax.grid(True)
    ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    #plt.yscale("log")
    #ax.tick_params(axis='y', labelrotation=90)
    ax.set_xlabel(X_varlabel)
    ax.set_ylabel(Y_varlabel)

    plt.savefig(f'{path_fig}{foldername}.{Yvarname}of{Xvarname}.{loc_cp}.png', dpi = 300)
    plt.close()

def plot_zofX(Xvarname, Xcp, tmod,**kwargs): # hysteretic curve

    #Xvarname = 'satur'

    X = Xcp[0] # single point
    loc_cp = Xcp[1]

    time_indexes = tidxscaling(dt_real,tmod,c=kwargs.get('c',float),d=kwargs.get('d',float),e=kwargs.get('e',float))
    #nt = len(time_indexes)
    #tarray = dt_real[time_indexes]

    # X
    if Xvarname == 'H':
        X_array_3Dt = bernoullihead_3Dt()
    else :
        vname = Xvarname
        X_array_3Dt = readpfblist_to_3Dtarray(path,runname,vname)
    X_array_Xt = X_array_3Dt[:,:,X[0],X[1]]
    pltsettings = variablescaling(X_array_3Dt, vname)
    X_varlabel = pltsettings[1]

    fig, ax = plt.subplots(1,2,figsize=(5,3.5), gridspec_kw={'width_ratios': [0.5,1],'height_ratios': [1]})
    
    ax[0].plot(permz[3:,X[0],X[1]],z_centers[3:])
    ax[0].set_xscale("log")
    ax[0].grid(True)
    ax[0].set_xlabel(r'K$_{sat}$ [m/h]')
    ax[0].set_xticks([1.0e-3,1.0e-2])
    ax[0].set_ylabel('z [m agl]')

    for idx in time_indexes:
        ax[1].plot(X_array_Xt[idx,3:],z_centers[3:],label=f'{round(dt_real[idx],0)}h', marker='.')
    
    ax[1].grid(True)
    ax[1].legend() #loc='upper right', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    #plt.yscale("log")
    #ax.tick_params(axis='y', labelrotation=90)
    ax[1].set_xlabel(X_varlabel)
    ax[1].set_ylabel('z [m agl]')

    plt.savefig(f'{path_fig}{foldername}.zof{Xvarname}.{loc_cp}z3-10.png', dpi = 300)
    plt.close()

################################################################
# Stability

# evaluation only at dump times
def CFL_eval():
    
    dt = Deltat # len = nt (number of dump times)

    #v = vel = velocitynorm_3Dt()
    vname = 'vx_c'
    data_3Dt = centered_var_generator(vname)
    vx = data_3Dt
    vname = 'vy_c'
    data_3Dt = centered_var_generator(vname)
    vy = data_3Dt
    vname = 'vz_c'
    data_3Dt = centered_var_generator(vname)
    vz = data_3Dt

    # cell centered variable
    CFL_norm_matrix = np.empty((nt,nz,ny,nx))
    CFL_x_matrix = np.empty((nt,nz,ny,nx))
    CFL_y_matrix = np.empty((nt,nz,ny,nx))
    CFL_z_matrix = np.empty((nt,nz,ny,nx))

    for t in range(nt) :
        #CFL_norm_matrix[t] = v[t]*dt[t]/dx_norm
        CFL_x_matrix[t] = vx[t]*dt[t]/dx
        CFL_y_matrix[t] = vy[t]*dt[t]/dy
        CFL_z_matrix[t] = vz[t]*dt[t]/dz_3D

    #CFL_norm_stat = [mean(), var(), max(CFL_norm_matrix),np.where(CFL_matrix = max(CFL_matrix.max))]
    #CFL_x_stat = [mean(), var(), max(CFL_norm_matrix),np.where(CFL_matrix = max(CFL_matrix.max))]
    #CFL_y_stat = [mean(), var(), max(CFL_norm_matrix),np.where(CFL_matrix = max(CFL_matrix.max))]
    #CFL_z_stat = [mean(), var(), max(CFL_norm_matrix),np.where(CFL_matrix = max(CFL_matrix.max))]

    #return CFL_x_stat, CFL_y_stat, CFL_z_stat
    return CFL_x_matrix, CFL_y_matrix, CFL_z_matrix

# max, mean, min ratios of dt/dX
dt_stat = [max(Deltat), np.mean(Deltat), min(Deltat)]
dz_max = max(dz)
dX_stat = [dz_max, dx, min(dz)]
ratio = [dt_stat[0]/dX_stat[2],dt_stat[1]/dX_stat[1],dt_stat[2]/dX_stat[1]]
print("ratio", ratio)

def plot_CFL_evol():

    CFL_X = CFL_eval()

    CFL_xt = maxinspace_3Dt(CFL_X[0])
    CFL_yt = maxinspace_3Dt(CFL_X[1])
    CFL_zt = maxinspace_3Dt(CFL_X[2])

    fig, ax = plt.subplots(1)

    ax.plot(dt_real, CFL_xt , label = f'CFL_x', marker='.', color='b')
    ax.plot(dt_real, CFL_yt, label = f'CFL_y', marker='.', color='r')
    ax.plot(dt_real, CFL_zt, label = f'CFL_z', marker='.', color='g')
    
    ax.grid(True)
    ax.legend(loc='upper right') #, bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    ax.set_xlabel('t [h]')
    ax.set_ylabel('CFL')

    plt.savefig(f'{path_fig}{foldername}.CFL.png', dpi = 300)

def plot_timestep():
    
    #tzoom = 70
    tzoom = -1

    fig, ax = plt.subplots(1,figsize=(5,2))
    ax.plot(dt_real[:tzoom], Deltat[:tzoom] , marker='.', color='k')
    
    ax.grid(True)
    #ax.legend(loc='upper right') #, bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    ax.set_xlabel('t [h]')
    ax.set_ylabel('timestep [h]')

    plt.savefig(f'{path_fig}{foldername}.dtstep.png', dpi = 300)

plot_timestep()

def val_derivativecv():

    vname = 'water_table'
    array_2Dt = -1* readpfblist_to_2Dtarray(path,runname,vname)

    ## FIND TIME OF CONVERGENCE
    # temporal derivative of wt
    deriv_wt = derivative2D_dvarofdt(array_2Dt)
    # mean in space of temporal derivative of wt
    mean_deriv_wt = layer_mean1Dt(deriv_wt)
    abs_mean_deriv_wt = abs(mean_deriv_wt)
    threshold = 5.0e-8
    idx = np.where(abs_mean_deriv_wt<threshold)
    idx = idx[0][0]
    
    ## ASSESS CONVERGENCE PROXY AT TIME OF CONVERGENCE
    array_known = unconfined_PBC_2D()
    diffres = diff_array2D(array_2Dt[idx],array_known)
    diff_2D = diffres[0]
    diff_1D = diff_2D[0,:,int(nx/2)]
    diffratio_1D = diff_1D/array_known[0,:,int(nx/2)]

    val_cv = abs(diffratio_1D.mean())
    
    return val_cv

################################################################
# THEORETICAL EQUATION

def unconfined_Dupuit(Q,R,K,ho,y_array): # steady state
    return (Q/(K*np.pi)*np.log10(y_array/R) + ho**2)**(1/2) + z_faces[0]

def unconfined_PBC(hup,hlow,Ly,y): #steady state
    return ((hup**2-hlow**2)/Ly*y+hlow**2)**(1/2) + z_faces[0]

def unconfined_PBC_2D(): #steady state
    array_2D = np.ones((1,ny,nx))

    hup = -1*z_faces[0] + press_BC['y-upper__alltime']
    hlow = -1*z_faces[0] + press_BC['y-lower__alltime']
    #print(hlow)
    array_Dupuit = unconfined_PBC(hup,hlow,ny*dy,y_centers)
    for i in range(nx):
        array_2D[0,:,i] *= array_Dupuit
    # array_2D[0,:] = array_2D[0,:]*array_Dupuit
    return array_2D

################################################################

# Nomenclature of variables

pfb_out3dtvariables = ['press','satur','velx','vely','velz']
pfb_cellcentered3dtvariables = ['press', 'satur', 'vx_c','vy_c','vz_c','v', 'H','h'] #'subsurface_storage' - H for hydraulic head
#pfb_inparams = [slopex, slopey, permx, permy, permz, manning, porosity, ... ] #nb slope can be 2d or 3d ?
pfb_out2dtvariables = ['water_table', 'overland_flow'] #'surface_storage', 'surface_storage', 'overland_bc_flux', 'overlandsum',
faces = ['x-lower', 'x-upper', 'y-lower', 'y-upper', 'z-lower', 'z-upper']
plot_projections = ['2dxy','2dxz','2dyz','1d','3d']


zidx = nz-1
yidx = 0
xidx = int(nx/2)

#plot_proj2d_multivardtall(runname,'2dxy',10, "int", c=7, d=dt_real[0], e=dt_real[8])


# plot_3d_geom()


tmod = "int" # all
#tmod = "beg" # beginning
#tmod = "int" # interval
#tmod = "ed"
ndt = 7
tinit = dt_real[0] # initial time in hours - for plot
tfin = 165 # dt_real[-1] # final time in hours - for plot
plot_proj2d_multivardtall(runname,'2dxz',yidx, tmod, c=ndt, d=tinit, e=tfin)
## plot_proj2d_multivardtall(runname,'2dxz',yidx, "all")

plot_proj2d_multivardtall(runname,'2dxy',zidx, tmod, c=ndt, d=tinit, e=tfin)
plot_proj2d_multivardtall(runname,'2dyz',xidx, tmod, c=ndt, d=tinit, e=tfin)


#plot_3d_singlevardtsgl('H', 0) # for cellcentered3dtvariables
#plot_3d_singlevardtsgl('H', 1) # for cellcentered3dtvariables
#plot_3d_singlevardtsgl('H', 2) # for cellcentered3dtvariables
#plot_3d_singlevardtsgl('H', 3) # for cellcentered3dtvariables
#plot_3d_singlevardtsgl('H', 30) # for cellcentered3dtvariables


#layerstoplot = [10,9,8,7,6,5,4,3,2,1,0]  #VM
layerstoplot = [5,4,3,2,1,0]
#plot_proj2d_singlevardtslall('press','2dxy',layerstoplot, tmod, c=ndt, d=tinit, e=tfin)
#plot_proj2d_singlevardtslall('velz','2dxy',layerstoplot, tmod, c=ndt, d=tinit, e=tfin)
#plot_proj2d_singlevardtslall('v','2dxy',layerstoplot, tmod, c=ndt, d=tinit, e=tfin)
#plot_proj2d_singlevardtslall('satur','2dxy',layerstoplot, tmod, c=ndt, d=tinit, e=tfin)
#plot_proj2d_singlevardtslall('h','2dxy',layerstoplot, tmod, c=ndt, d=tinit, e=tfin)

## VM check points
f_cp = '/home/patras/PF-Valmalenco/data/controlpoints.txt'
Xidx_cp = read_cpcsv(f_cp)
Xidx_cp0 = [[Xidx_cp[0][0,0],Xidx_cp[0][0,1]],Xidx_cp[1][0]]
Xidx_cp1 = [[Xidx_cp[0][1,0],Xidx_cp[0][1,1]],Xidx_cp[1][1]]
Xidx_cp2 = [[Xidx_cp[0][2,0],Xidx_cp[0][2,1]],Xidx_cp[1][2]]
#print('Control Points', Xidx_cp)
print('check p 0',Xidx_cp0[0])

print('slope in Palu', z_dem_2D[Xidx_cp0[0][0],Xidx_cp0[0][1]], slopex[0,Xidx_cp0[0][0],Xidx_cp0[0][1]], slopey[0,Xidx_cp0[0][0],Xidx_cp0[0][1]]) #,
print('slope in Gombaro', z_dem_2D[Xidx_cp1[0][0],Xidx_cp1[0][1]], slopex[0,Xidx_cp1[0][0],Xidx_cp1[0][1]], slopey[0,Xidx_cp1[0][0],Xidx_cp1[0][1]])
print('slope in Spriana', z_dem_2D[Xidx_cp2[0][0],Xidx_cp2[0][1]], slopex[0,Xidx_cp2[0][0],Xidx_cp2[0][1]], slopey[0,Xidx_cp2[0][0],Xidx_cp2[0][1]])

XP_ylower = [Xidx_cp[0][1],[Xidx_cp[1][1]]]
XP_center = [np.array([int(ny/2),int(nx/2)]),f'[{int(ny/2)},{int(nx/2)}]']

## 'Universal' check point
#Xidx_cp = [np.array([[0,int(nx/2)],[1,int(nx/2)],[2,int(nx/2)],[int(ny/2),int(nx/2)],[ny-1,int(nx/5)],[ny-1,int(nx/2)]]),['P1','P2','P3','P4','P5','P6']]
XP_center = [np.array([int(ny/2),int(nx/2)]),f'[{int(ny/2)},{int(nx/2)}]']
#XP_ylower = [np.array([0,int(nx/2)]),f'[0,{int(nx/2)}]']
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


tmod = "intall"
#tmod = "beg" # beginning
#tmod = "intall" # interval
ndt = 7
tinit = 0 # initial time in hours - for plot
tfin = dt_real[nt-1] # final time in hours - for plot
"""
plot_ZXoft(runname,Xidx_cp,tmod, c=ndt, d=tinit, e=tfin) # X : multiple points

plot_Zxoft(runname,Xidx_cp0,tmod, c=ndt, d=tinit, e=tfin)
plot_Zxoft(runname,Xidx_cp1,tmod, c=ndt, d=tinit, e=tfin)
plot_Zxoft(runname,Xidx_cp2,tmod, c=ndt, d=tinit, e=tfin)

plot_Zxoft(runname,XP_ylower,tmod, c=ndt, d=tinit, e=tfin) # x : single point
plot_Zxoft(runname,XP_center,tmod, c=ndt, d=tinit, e=tfin)

plot_Hxoft(runname,Xidx_cp0,tmod, c=ndt, d=tinit, e=tfin)
plot_Hxoft(runname,Xidx_cp1,tmod, c=ndt, d=tinit, e=tfin)
plot_Hxoft(runname,Xidx_cp2,tmod, c=ndt, d=tinit, e=tfin)

plot_Hxoft(runname,XP_ylower,tmod, c=ndt, d=tinit, e=tfin) # Hxoft = H(x)(t), evolution of H(x) in time
plot_Hxoft(runname,XP_ylowerbis,tmod, c=ndt, d=tinit, e=tfin)
plot_Hxoft(runname,XP_center,tmod, c=ndt, d=tinit, e=tfin)
"""
# plot_compared_surfaceflow(runname,XP_maxof,tmod, d=tinit, e=tfin)

plot_compared_surfaceflow(runname,Xidx_cp0,tmod, c=ndt, d=tinit, e=tfin)
plot_compared_surfaceflow(runname,Xidx_cp1,tmod, c=ndt, d=tinit, e=tfin)
plot_compared_surfaceflow(runname,Xidx_cp2,tmod, c=ndt, d=tinit, e=tfin)

plot_boundarychecks(tmod,c=ndt, d=tinit, e=tfin)

plot_YofX('h','satur', Xidx_cp0)
plot_YofX('h','satur', Xidx_cp1)
plot_YofX('h','satur', Xidx_cp2)

tmod = "ed"
#tmod = "beg" # beginning
#tmod = "intall" # interval
ndt = 7
tinit = 0 # initial time in hours - for plot
tfin = dt_real[nt-1]
plot_zofX('satur', Xidx_cp2, tmod, c=ndt, d=tinit, e=tfin)
plot_zofX('satur', Xidx_cp1, tmod, c=ndt, d=tinit, e=tfin)
plot_zofX('satur', Xidx_cp0, tmod, c=ndt, d=tinit, e=tfin)
plot_zofX('press', Xidx_cp2, tmod, c=ndt, d=tinit, e=tfin)

#tmod = "fin"
ndt = 7
tmod = "int" # interval
tinit = dt_real[0] # initial time in hours - for plot
tfin = dt_real[-1]-1 # 5e7 # final time in hours - for plot #print(dt_real[-1])
#plot_HTofy(runname,tmod, c=ndt, d=tinit, e=tfin ) # H(t) distributed in space along y axis (x=nx/2)

cv_ratio = 0.1
derive_threshold = 5e-7

#plot_MXFig3()

if runname == 'DSc1000z10s0':
    #plot_HTofybis()
    ndt = 6
    tmod = "int" # interval
    tinit = dt_real[0] # initial time in hours - for plot
    tfin = dt_real[-1]-1 #6000 #
    plot_DHTofy_toknownsolution(tmod, c=ndt, d=tinit, e=tfin)
    tmod = 'all'
    cvproxy_toknownsolution(tmod,cv_ratio)

plot_CFL_evol()

# v46=v45_ time of cv (<5%) t[10] = 6311.923
# v42 time of cv (<5%) t[13] = 10911980.0

plot_wtderivative()

#val = val_derivativecv()
#print('RelativeCV proxy at t|dh/dt<5e-8:', val)
