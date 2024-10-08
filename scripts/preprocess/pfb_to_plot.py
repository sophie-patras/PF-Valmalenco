#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
Pfb file (2d (1-layer) or 3d (p-layers)) to stacked plot

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""

# import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb

parser = argparse.ArgumentParser()
parser.add_argument('pfbfilename', type = str)
args = parser.parse_args()

fig_dir = "/mnt/c/Users/Sophie/Documents/4-Figures/"

#path_out = "/home/patras/Valmalenco/Data/DataPF/"
#filename_out = "slopeY.c500.v2.pfb"
# path_out = "/home/patras/Lombardy/Tmp/PLT_35_Leo/" #pyyamlshTmp/"
# path_out = "/home/patras/Valmalenco/Tmp/PLT/"
# path_out = "~/Sabino_grapp1-mb_sensitivity-7835310/model_runs/PFCLM/"
path = "../../data/pfdata/"
# filename_out = "PLT_Box.out.press.00000.pfb"
# filename_out = "PLT_Box.out.slope_x.pfb"
filename_out = args.pfbfilename 
f_out = path + filename_out
data_to_plot = read_pfb(get_absolute_path(f_out))

datashape = data_to_plot.shape
print(f'(NZ,NY,NX): {datashape}') # plot (NZ,NY,NX)
#print(type(data_to_plot))  # numpy array

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))
#y = y.max() - y

print('NB. Only 3 layers ploted')
dl = int(datashape[0]/4)+1

for p in range(datashape[0]):
    data = data_to_plot[p]
    dmin = data.min()
    dmax = data.max()
    print('Layer',p,'- min,max:',dmin,',',dmax)
    if p in range(0,datashape[0],dl):
        # use mask to avoid e38 value. https://sigon.gitlab.io/post/2018-11-08-plot-raster-nodata-values/
        ax.plot_surface(x,y,np.full_like(data, p), facecolors=plt.cm.viridis(data), rstride=1, cstride=1, antialiased=True, shade=False)

ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_zlabel('i')
ax.set_title(filename_out[:-4])

m = plt.cm.ScalarMappable(cmap=plt.cm.viridis)  #, norm=surf.norm)
#dmin = data.min()
#dmax = data.max()
#ax.set_zlim(vmin,vmax)
#m.set_clim(max(-1,vmin),vmax)
m.set_clim(dmin,dmax)
#print('min,max:',dmin,',',dmax)
#print(f'non all nul in cell : {np.where(data>-3.402823466385288e+38)}')
plt.colorbar(m, ax=plt.gca())

plt.show()
plt.savefig(fig_dir+args.pfbfilename[:-4]+'.png')
#plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
