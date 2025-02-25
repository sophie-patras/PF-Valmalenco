#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For graphical output
Convert real crs coordinates, to computational closest cell centered coordinates (es. application for control points, in controlpoint.csv - EPSG:32632)
    - read asc dem (in EPSG:32632)
    - location of the control point in (i,j) coordinates

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""
# import os
import numpy as np
import pandas as pd

###############################################################################################
# INPUT
# input .asc dem
path_dem = '/home/patras/PF-Valmalenco/data/prepareddata/'
filename_dem = 'hydroDEM.c250.v2.asc'
f_dem = path_dem + filename_dem
# controlpoints
f_cp = '/home/patras/PF-Valmalenco/data/controlpoints.txt'
###############################################################################################

def read_demheader(f_dem):
    
    header_rows = 6
    header_info = {}
    row_ite = 1
    with open(f_dem, 'rt') as file_h:
        for line in file_h:
            if row_ite <= header_rows:
                line = line.split(" ", 1)
                header_info[line[0]] = float(line[1])
            else:
                break
            row_ite = row_ite+1
    
    # read data array :: dem
    # dem = np.loadtxt(f_dem, skiprows=header_rows, dtype='float64')
    # print(dem)

    # EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
    ncols = int(header_info['ncols'])
    nrows = int(header_info['nrows'])
    cellsize = header_info['cellsize']
    xwest = header_info['xllcorner'] # x lower
    xright = xwest+ncols*cellsize
    ysouth = header_info['yllcorner'] # y lower
    ynorth = ysouth+nrows*cellsize

    return ncols,nrows,cellsize,xwest,xright,ysouth,ynorth

def read_cpcsv(f_cp):
    # csv containing list of control points

    data_cp = pd.read_csv(f_cp)
    x_cp = np.array(data_cp['xcoord'])
    #print(x_cp)
    y_cp = np.array(data_cp['ycoord'])
    #print(y_cp)
    loc_cp = data_cp['localita']

    # pfb read from lower left
    demheader = read_demheader(f_dem)
    cellsize = demheader[2]
    ysouth = demheader[5]
    xwest = demheader[3]
    i_cp = [int((y_cp[k] - (ysouth-cellsize/2))/cellsize) for k in range(3)]
    j_cp = [int((x_cp[k] - (xwest-cellsize/2))/cellsize) for k in range(3)]

    CP = np.array([i_cp,j_cp])
    #print(CP)
    #print(CP.transpose())
    #Xidx = np.reshape([i_cp,j_cp],(len(i_cp),2))
    Xidx = CP.transpose()

    return Xidx, loc_cp

#Xidx = read_cpcsv(f_cp)
#print(Xidx[0].shape)
#print(Xidx)