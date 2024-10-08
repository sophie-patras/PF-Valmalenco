#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
For file format conversion, from .asc (6 lines header, first data: north west) to .sa - used in parflow - (1 line header, first data: south west)

Created on Sat Sep 7 18:02:20 2024
@author: S.P.
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ascfilename', type = str)
args = parser.parse_args()

# input .asc
path_in = '/home/patras/Valmalenco/Data/DataElab/'
#path_in = '/home/patras/Lombardy/Data/TXT/'
filename_in = args.ascfilename
#filename_in = 'hydroDEM.c500.v2.asc'
f_in = path_in + filename_in
f_out = path_in + filename_in[:-4] + '.sa'

# read
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
# read data array :: dem
dem = np.loadtxt(f_in, skiprows=header_rows, dtype='float64')
# print(dem)

# EPSG: 32636 :: WGS84UTM32N - cartesian - unit [m]
ncols = int(header_info['ncols'])
nrows = int(header_info['nrows'])

# write
fout = open(f_out,"w")
fout.write(str(ncols)+" "+str(nrows)+" 1\n")
for i in range(nrows):
    for j in range(ncols):
        fout.write(str(dem[nrows-1-i,j])+"\n")
fout.close()
