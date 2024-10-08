# -*- coding: utf-8 -*-

"""
Read .tif file in georeferenced coordinate system

Created on Mon Oct 7 10:47:18 2024
@author: SP.
"""

#import numpy as np
import matplotlib.pyplot as plt
#I = plt.imread(tiff_file)

import rioxarray
import rasterio

path_in = "../../data/rawdata/SG_S14_2017_ESDAC/"
fname_in = "KS_M_sl4_EU_047_014.tif"
f_tif = path_in + fname_in

path_out = "../../data/prepareddata/"
#fname_out = fname_in[:-4] + ".32632" + ".tif"

fig_dir = "/mnt/c/Users/Sophie/Documents/4-Figures/"

band = rioxarray.open_rasterio(f_tif)
print(type(band))
#band.plot()
#plt.savefig(fig_dir + fname[:-4] + ".77261.png")
print(band.rio.crs,' original')

#Reproject the raster to geographic coordinates
band_geographic = band.rio.reproject('EPSG:32632')
print(band_geographic.rio.crs)
#band_geographic.plot()
#plt.savefig(fig_dir + fname[:-4] + ".32632.png")

#Reproject the raster to UTM coordinates
#band_utm = band_geographic.rio.reproject(band_geographic.rio.estimate_utm_crs())

#band_geographic.rio.to_raster(path_out) not working
