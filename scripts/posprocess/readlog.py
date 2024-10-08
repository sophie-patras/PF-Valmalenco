#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Read <runname>.out.log files
    - read to array {'Sequence','Time','\Delta t','Dumpfile'}
    - 

Created on Tue Sep 03 13:34:20 2024
@author: S.P.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

###
np.set_printoptions(precision=6)
FIG_DIR = '/mnt/c/Users/Sophie/Documents/4-Figures/'
plt.style.use('~/PF-Valmalenco/Codes/settings.mplstyle')
###


#### READ PFB output results

parser = argparse.ArgumentParser()
parser.add_argument('runname', type = str)
args = parser.parse_args()

# Input log
path_in = './'
fn_in = args.runname+'.out.log'
f_in = path_in + fn_in

# Lines to read:
#"Total Timesteps: 1277" (l12)
#"Total Run Time: 9659.799000 seconds" (l2660)

skip_rows1 = 11 #then read tot timesteps
skip_rows2 = 4 #loadtxt

header_rows = skip_rows1 + 1 + skip_rows2
header_info = {}
row_ite = 1
with open(f_in, 'rt') as file_h:
     for line in file_h:
        if row_ite == 12:
             line = line.split(":", 1)
             header_info[line[0]] = int(line[1])
        elif row_ite > 12:
             break
        row_ite = row_ite+1
# read data array :: dem

sequencesnumber = header_info['Total Timesteps']
#print('Total TimeSteps: ', sequencesnumber)
print(header_info)

dataset = np.loadtxt(f_in, \
        skiprows=header_rows, \
        max_rows= sequencesnumber, \
        usecols=(0, 1, 2, 4), \
        dtype={'names': ('sequence', 'time', 'timestep', 'dumptime'), \
        'formats': ('float64', 'float64', 'float64', '|S15')})

time = dataset['time']
dumptime = np.array(dataset['dumptime'])
#print(dumptime)
#print(len(dumptime))
idxdt = np.where(dumptime!=b'y')[0]
#print(idxdt)
time = time[idxdt]
dumptime = np.linspace(0,len(idxdt),len(idxdt))

#for s in range(sequencesnumber):
#    if str(dumptime[s])=="b'y'":
#        dumptime[s]= None #np.nan

fig = plt.figure(1,figsize=(4.45,4.45))
ax = fig.add_subplot(111)

ax.plot(time,dumptime,marker='+',linestyle='')
ax.grid(True)
ax.set_xlabel('t [h]')
ax.set_ylabel('dumpfile nb')

plt.savefig(FIG_DIR + args.runname + '.log.png')


def dumprealtime(dumpnb):
    return time[dumpnb]

print("real time of dumpfile '00007' : ",dumprealtime(7),"h")
#print("real time of dumpfile '00030' : ",dumprealtime(30),"h")
