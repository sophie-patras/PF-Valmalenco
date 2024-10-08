# -*- coding: utf-8 -*-
"""
PreProcessing
Times series from meteorological monitoring stations

Created on Mon Mar 25 16:28:15 2024
@author: SP.
"""
import numpy as np
import pandas as pd
import re
import os
import argparse

import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

marker_style = {'marker':'.','markersize':'8'}
line_style = {'linestyle':'-','linewidth':'1.5'}
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.family'] =['URW Bookman'] #Perpetua

# https://earthly.dev/blog/plotting-rainfall-data-with-python-and-matplotlib/
Noval = -9999 #111  #No value (usually -9999)

parser = argparse.ArgumentParser()
parser.add_argument('tipology',type = str) # pluvio, temp, ur, windir, windvel, rad, snow
#parser.add_argument('dumpt',type = int)
#parser.add_argument('valmaxplot',type = int)
# parser.add_argument('tipology', type = str) # pluvio, temp, ur, windir, windvel, rad, snow
args = parser.parse_args()

# DB Keys ref from 'AnagraficaSensoriWeb.csv'

Tipologia = [['Precipitazioni',"mm"],\
        ["Temperatura","°C"],\
        ["Umidità Relativa","%"],\
        ['Direzione vento',"°"],\
        ["Velocità vento","m/s"],\
        ["Radiazione Globale","W/m^2"],\
        ["Altezza neve","cm"],\
        ["Livello idrometrico","cm"]]

NomeStazione = ["Alpe Entova (Chiesa in Valmalenco)",\
        "Passo Marinelli (Lanzada)",\
        "Palù (Lanzada)",\
        "Ganda (Lanzada)",\
        "Piazzo Cavalli (Caspoggio)",\
        "Fond. Fojanini (Sondrio)"]

IdSensore = np.array([[41,53,12754,Noval,Noval,Noval,12753,Noval],\
              [14224,14221,14218,14216,14217,14207,14205,Noval],\
              [8094,8092,32259,9763,11754,32363,8093,Noval],\
              [46,56,17545,Noval,Noval,Noval,Noval,64],\
              [48,58,10773,132,131,Noval,70,Noval],\
              [2103,2096,2097,2104,11736,2098,Noval,Noval]])

def yyyymmdd_parser(s):
    return datetime.datetime.strptime(s, '%Y/%m/%d %H:%M')


# Merging multiple datasets
def get_data(weather_dir):
    # p = re.compile('data_[1]+.csv')
    p = re.compile('RW_'+idsensore+'_20[1-2][0-9]+.csv')
    data = []
    for filename in sorted(os.listdir(weather_dir)):
        if p.match(filename):
            path = os.path.join(weather_dir, filename)
            # We identify which column contains the date, along with the date parser to use.
            yearly_data = pd.read_csv(path, parse_dates=['Data-Ora'], date_parser=yyyymmdd_parser)
            data.append(yearly_data)
    return pd.concat(data)

#function that splits up the data by year, and calculates a new field that is days since January 1st of the current year

def get_days_since_jan1(date):
    return date.timetuple().tm_yday - 1

def split_data_by_year(weather):
    yearly_rainfall = []
    yearly_days_since_jan1 = []
    current_year = None
    by_year = []
    for i, row in weather.iterrows():
        date = row['Date/Time']
        rainfall = row['Total Precip (mm)']

        if current_year is None:
            current_year = date.year
        elif current_year != date.year:
            by_year.append((current_year, pandas.DataFrame({
                'days_since_jan1': yearly_days_since_jan1,
                'rainfall': yearly_rainfall,
                })))
            yearly_days_since_jan1 = []
            yearly_rainfall = []
            current_year = date.year

        yearly_rainfall.append(rainfall)
        yearly_days_since_jan1.append(get_days_since_jan1(date))

    if current_year is not None:
        by_year.append((current_year, pandas.DataFrame({
            'days_since_jan1': yearly_days_since_jan1,
            'rainfall': yearly_rainfall,
            })))

    return by_year

###################################################################################

tipo = args.tipology
datadir = "/home/patras/Valmalenco/Data/DataRaw/METEO_VM_20XX_AL/"
dumptime = 60 #args.dumpt # in {10;30; 60; 1440} in min


if (tipo=='pluvio'):

    plotnb = 6 #len(IdSensore[:,0])
    fig, ax = plt.subplots(2,3,figsize=(15.0,10.0))
    ploti = 1
    
    for idsensore in IdSensore[:,0]:
        idsensore = str(idsensore)
        dataset = get_data(datadir)

        values = np.array(dataset['Valore Cumulato'])
        date = [x.to_pydatetime() for x in dataset['Data-Ora']]

        #idex_ok = np.where(values >= -1)[0]
        idex_nan = np.where(values <= -1)[0]
        #print('Warning : no data values \n', idex_nan)
        values[idex_nan] = 0 #or use pandas fillna(0)
        #values = np.array(values)[idex_ok]
        #date = np.array(date)[idex_ok]

        # correct 'border' no data
        values[idex_nan]=(values[idex_nan-1]+values[idex_nan+1])/2
        print('For /10min measure')
        print('mean', np.mean(values))
        print('min', values.min())
        print('max', values.max())

        if dumptime>10:
            dumpint = int(dumptime/10)
            nbval = int(len(date)/dumpint)
            dumpdate = [date[i] for i in range (dumpint,len(date),dumpint)]
            dumpvalues = [sum(values[dumpint*i:dumpint*(i+1)-1]) for i in range (0,nbval)]
            date = dumpdate
            values = dumpvalues

            print('For /60min measure')
            print('mean', np.mean(values))
            print('min', min(values))
            print('max', max(values))

        ax = fig.add_subplot(2,3,ploti)        
        ax.plot(date,values)
        ax.set_ylim(ymin=0, ymax=35) #args.valmaxplot) # 35 for 1h
        ax.set_yticks(np.linspace(0,35,8))
        #ax.set_xlabel('date')
        # ax.set_ylabel('precipitation (mm)')
        ax.set_title('sensor '+ idsensore)
        ax.grid(True)
        ploti = ploti + 1

#plt.title('Cumulated precipitation in mm/'+ str(dumptime) + 'min')
plt.show()

#plotdir = ""
#plotname = 'chart_'+idsensore+'_per'+str(dumptime)+'min.png'
#fig.savefig(plotdir+plotname, bbox_inches='tight')
