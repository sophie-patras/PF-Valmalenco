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
# np.set_printoptions(precision=6)

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant
###

###############################################################################################
# INPUT

## MX
path = '/home/patras/PF-Test/Maxwell2013/outputs/'
foldername = 'MX.c1s1y3_v6' #'MX.c100s1bcx_v5'
runname = 'MX.c1s1'

## DS
path = '/home/patras/PF-Test/DumbSquare/outputs/'
foldername = 'DSc1000z10s0.MX.DP23.IN0dt01_v55' #'DS.c100s1_v28' DSc10z10s0.MX.DP23.IN0_v52
runname = 'DSc1000z10s0'

## VM
#path = '/home/patras/PF-Valmalenco/outputs/'
#foldername = 'CLM_V52'
#runname = 'CLM_V5'

## VM
#path = '/home/patras/PF-Valmalenco/outputs/'
#foldername = 'SeepBox_KSsgL11IC-0125' #'CLM_V52'
#runname = 'SeepBox'

#path = '/home/patras/PF-Valmalenco/outputs/BowlBox_UXZ/'
#foldername = 'BowlBox_UXZ_K36e-3L6d' #'CLM_V52'
#runname = 'BowlBox_UXZ'

## VM.RR
path = '/home/patras/PF-Valmalenco/outputs/'
foldername = 'RainRec_RFe-2C2-2IC-35'
runname = 'RainRec_RFe-2C2-2IC-35'

## LW
#path = '/home/patras/PF-Test/LW/outputs/'
#foldername = 'LW_var_dz_spinup'
#runname = 'LW_var_dz_spinup'

###############################################################################################

# Input log

def log_to_array(path, foldername, runname):
     
     f_in = f'{path}{foldername}/{runname}.out.log'

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

     sequencesnumber = header_info['Total Timesteps']
     #print(header_info)

     dataset = np.loadtxt(f_in, \
          skiprows=header_rows, \
          max_rows= sequencesnumber+1, \
          usecols=(0, 1, 2, 4), \
          dtype={'names': ('sequence', 'time', 'timestep', 'dumptime'), \
          'formats': ('float64', 'float64', 'float64', '|S15')})
     
     # turn to no data the non dumpedtimes

     return dataset

def dump_to_simulatedtimes_equivalent(path, foldername, runname):

     dataset = log_to_array(path, foldername, runname)
     #print(dataset)

     time = np.array(dataset['time'])
     dumptime = np.array(dataset['dumptime'])
     idxdt = np.where(dumptime!=b'y')[0]
     #print(idxdt)

     #dumptimes = np.linspace(0,len(idxdt)-1,len(idxdt))
     #print(dumptimes)
     #dumpdataset = np.array([dumptimes, time[idxdt]])
     if len(idxdt)>0:
          dumpdataset = time[idxdt]
     else:
          dumpdataset = time

     return dumpdataset

#dt_real = dump_to_simulatedtimes_equivalent(path, foldername, runname)
#print(dt_real)

def dump_timesteps(path, foldername, runname):

     dataset = log_to_array(path, foldername, runname)
     #print(dataset)

     timestep = np.array(dataset['timestep'])
     dumptime = np.array(dataset['dumptime'])
     idxdt = np.where(dumptime!=b'y')[0]
     #print(idxdt)

     #dumptimes = np.linspace(0,len(idxdt)-1,len(idxdt))
     #print(dumptimes)
     #dumpdataset = np.array([dumptimes, time[idxdt]])
     if len(idxdt)>0: # solution if timesteps gave outputs pfb
          dumpdataset = timestep[idxdt]
     else:
          dumpdataset = timestep

     return dumpdataset

"""
def nt_dump():
    return len(idxdt)+1

print(nt_dump())

def dt_real():
     return time[idxdt]

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

plt.savefig(f'{path_fig}{foldername}.log.png')

def dumprealtime(dumpnb):
    return time[dumpnb]

print("real time of dumpfile '00007' : ",dumprealtime(7),"h")
#print("real time of dumpfile '00030' : ",dumprealtime(30),"h")
"""

dataset = log_to_array(path, foldername, runname)
sequence = np.array(dataset['sequence'])
timestep = np.array(dataset['timestep'])
time = np.array(dataset['time'])

timesteps_dumped = dump_timesteps(path, foldername, runname)
time_dumped = dump_to_simulatedtimes_equivalent(path, foldername, runname)

fig = plt.figure(1,figsize=(4.45,4.45))
ax = fig.add_subplot(111)

ax.plot(sequence,timestep,marker='+',linestyle='')
ax.grid(True)
ax.set_xlabel('sequence [#]')
ax.set_ylabel('dt [h]')

plt.savefig(f'{path_fig}{foldername}.dtsqc.log.png')

###
fig = plt.figure(2,figsize=(4.45,4.45))
ax = fig.add_subplot(111)

ax.plot(time,timestep,marker='+',linestyle='')
ax.plot(time_dumped,timesteps_dumped,marker='+', color='k', linestyle='',label='dumptimes')
ax.grid(True)
ax.legend()
ax.set_xlabel('t [h]')
ax.set_ylabel('dt [h]')

plt.savefig(f'{path_fig}{foldername}.dttime.log.png')
