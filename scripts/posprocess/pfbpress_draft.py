"""
For PF postprocessing
Pressure evolution

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""

# import os
import numpy as np
import matplotlib.pyplot as plt

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

# plotParams
# text_kwargs = dict(ha='center', va='center', fontsize=12, color='k')

# Spatial extent
# input DEM
path_in = '/home/patras/Valmalenco/Data/DataElab/'
filename_in = 'hydroDEMc500v2EPSG32632.asc'
f_in = path_in + filename_in

header_rows = 6
header_info = {}
row_ite = 1
with open(f_in, 'rt') as file_h:
     for line in file_h:
        if row_ite <= header_rows:
             line = line.split(" ", 1)
             header_info[line[0]] = float(line[1])
        else:
             break
        row_ite = row_ite+1

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])

cellsize = header_info['cellsize']
xwest = header_info['xllcorner']
xright = header_info['xllcorner']+ncols*cellsize
ysouth = header_info['yllcorner']
ynorth = header_info['yllcorner']+nrows*cellsize

cellsize = 2000
xwest = 0
ysouth = 0
# map_extent = (xwest, xright, ysouth, ynorth)

# Temp extent
# dumpinterval from: 00000 to 00073
DumpGap = 120 #hours (ie. 5 days)

simudir = "PLT_new"

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

for DumpInt in range(3,4):

    fn = "/home/patras/Valmalenco/Tmp/"+simudir+"/PLT.out.press.0000" + str(DumpInt) + ".pfb"
    npdata = read_pfb(get_absolute_path(fn))
    #npdata = npdata[:,:,2:]

    # print(f'Dimensions of output file: {press_data.shape}') # plot (NZ,NY,NX)
    # print(type(press_data))  # numpy array
    
    datashape = npdata.shape

    x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
    x = xwest + x*cellsize
    y = ysouth + y*cellsize

    #ax.clear()

    for p in range(datashape[0]):
	    data = npdata[p]
	    ax.plot_surface(x,y,np.full_like(data, p), facecolors=plt.cm.viridis(data), rstride=1, cstride=1, antialiased=True, shade=False)

    print(data)

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('p')
    ax.set_title('stacked hydraulic head')

    m = plt.cm.ScalarMappable(cmap=plt.cm.viridis)  #, norm=surf.norm)
    vmin = data.min()
    vmean = data.mean()
    vmax = data.max()
    print(vmin,vmean,vmax)
    m.set_clim(vmin,vmax)
    # m.set_clim (vmean-2,vmean+2)

    plt.colorbar(m, ax=plt.gca())
    ax.text2D(0.1,0.9,f"time={DumpInt*DumpGap}h",transform = ax.transAxes)
    plt.show()
    #plt.pause(0.5)

# Check differences:

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')

for DumpInt in range(3,4):

    fn = "/home/patras/Valmalenco/Tmp/"+simudir+"/PLT.out.press.0000" + str(30) + ".pfb"
    npdata = read_pfb(get_absolute_path(fn))
    #npdata = npdata[:,:,2:]

    fn1 = "/home/patras/Valmalenco/Tmp/"+simudir+"/PLT.out.press.0000" + str(31) + ".pfb"
    npdata1 = read_pfb(get_absolute_path(fn1))
    #npdata1 = npdata1[:,:,2:]

    # print(f'Dimensions of output file: {press_data.shape}') # plot (NZ,NY,NX)
    # print(type(press_data))  # numpy array

    datashape = npdata.shape

    x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
    x = xwest + x*cellsize
    y = ysouth + y*cellsize

    #ax.clear()                                                                                                                      
    for p in range(datashape[0]):
        data = (npdata1[p]-npdata[p])
        ax.plot_surface(x,y,np.full_like(data, p), facecolors=plt.cm.viridis(data), rstride=1, cstride=1, antialiased=True, shade=False)
    print(data)

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('p')
    ax.set_title('stacked hydraulic head')

    m = plt.cm.ScalarMappable(cmap=plt.cm.viridis)  #, norm=surf.norm)
    vmin = data.min()
    vmax = data.max()
    m.set_clim(vmin,vmax)

    plt.colorbar(m, ax=plt.gca())
    ax.text2D(0.1,0.9,f"diff dt={1*DumpGap}h",transform = ax.transAxes)
    plt.show()


# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
