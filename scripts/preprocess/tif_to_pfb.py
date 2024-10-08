# -*- coding: utf-8 -*-

"""
Read .tif file in georeferenced coordinate system

Created on Mon Oct 7 10:47:18 2024
@author: SP.
"""

import re, os, argparse

import matplotlib.pyplot as plt
import numpy as np
import rioxarray
import rasterio
from shapely.geometry import box
import geopandas as gpd

from parflow.tools.io import write_pfb

plt.style.use('../config.mplstyle')


# Inputs %%%%%%%%%%%%%%%%%%%%%%%
path_tif = "../../data/rawdata/SG_S14_2017_ESDAC/"
param = "KS" # in {KS, THS, MRC}
#fname_tif = "_M_sl[0-7]_EU_047_014.tif" #4
#f_tif = path_tif + var + fname_tif

epsg = 'EPSG:32632'
#f_mask = "../../data/prepareddata/mask."
# using printed output of ascdem extent
# map_extent = (xwest, xeast, ysouth, ynorth)
# map_extent = (553915.4983940837, 578415.4983940837, 5112351.329324091, 5137851.329324091) 

path_out = "../../data/pfdata/"
fname_out = param + "c250.pfb"
f_out = path_out + fname_out

fig_dir = "/mnt/c/Users/Sophie/Documents/4-Figures/"

cellsize = 250.0
extent = [553915.4983940837, 5112351.329324091, 578415.4983940837, 5137851.329324091]
nlayers = 1
nrows = 98
ncols = 102

# Processing %%%%%%%%%%%%%%%%%%%%%%%
bbox = box(*extent)
geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=epsg)

# Merging multiple datasets
def get_data(tif_dir):
    p = re.compile(param + '_M_sl[0-7]_EU_047_014.tif')
    data = [] #np.empty((nlayers,ncols,nrows))
    for filename in sorted(os.listdir(tif_dir)):
        if p.match(filename):
            f_tif = os.path.join(tif_dir, filename)
            band = rioxarray.open_rasterio(f_tif)
            print(band.rio.crs,' original')
            if band.rio.crs != epsg:
                band = band.rio.reproject(epsg)
            clipped = band.rio.clip(geo.geometry, geo.crs)
            data.append(clipped) #[-l] = clipped
    return data

data = get_data(path_tif)
print(data)

datasurf = data[-1]
datasurf.plot()
plt.savefig(fig_dir + "KS.32632.png")
# clip to mapextent

#data = np.flip(np.transpose(data.reshape(nlayers, ncols, nrows), (0, 2, 1)), axis = 0)
datanp = np.float64(data[0])
plt.savefig(fig_dir + ".32632.png")

datanp = np.transpose(datanp,(0, 2, 1))
# write to stacked pfb
write_pfb(f_out, datanp, p=2, q=2, r=1, x=0, y=0, z=0, dx=cellsize, dy=cellsize, dz=1.0, z_first=True, dist=False)
