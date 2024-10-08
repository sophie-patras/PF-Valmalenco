# -*- coding: utf-8 -*-
"""
PreProcessing
Times series from meteorological monitoring stations

Created on Mon Mar 25 16:28:15 2024
@author: SP.
"""
import numpy as np
import pandas as pd
import re, os, argparse

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
parser.add_argument('idsensor',type = int) # in IdSensore
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
        "Fond. Fojanini (Sondrio)", \
        "Ponte Eifell (Sondrio)"]

IdSensore = np.array([[41,53,12754,Noval,Noval,Noval,12753,Noval],\
              [14224,14221,14218,14216,14217,14207,14205,Noval],\
              [8094,8092,32259,9763,11754,32363,8093,Noval],\
              [46,56,17545,Noval,Noval,Noval,Noval,64],\
              [48,58,10773,132,131,Noval,70,Noval],\
              [2103,2096,2097,2104,11736,2098,Noval,Noval],\
              [Noval,Noval,Noval,Noval,Noval,Noval,Noval,16]])

# WARNING : files do not include 31/12 data

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

def get_days_since_jan1(date):
    return date.timetuple().tm_yday - 1

def split_data_by_year(weather):
    yearly_rainfall = []
    yearly_days_since_jan1 = []
    current_year = None
    by_year = []
    for i, row in weather.iterrows():
        date = row['Data-Ora']
        rainfall = row['Valore Cumulato']

        if current_year is None:
            current_year = date.year
        elif current_year != date.year:
            by_year.append((current_year, pd.DataFrame({
                'days_since_jan1': yearly_days_since_jan1,
                'rainfall': yearly_rainfall,
                })))
            yearly_days_since_jan1 = []
            yearly_rainfall = []
            current_year = date.year

        yearly_rainfall.append(rainfall)
        yearly_days_since_jan1.append(get_days_since_jan1(date))

    if current_year is not None:
        by_year.append((current_year, pd.DataFrame({
            'days_since_jan1': yearly_days_since_jan1,
            'rainfall': yearly_rainfall,
            })))

    return by_year

#######################################################################################"

idsensoreint = args.idsensor

idsensorecheck = len(np.where(IdSensore==idsensoreint)[0])
if (idsensorecheck==0):
    print('Error: idsensor not in list')
    quit()

idsensore = str(idsensoreint)
datadir = "/home/patras/Valmalenco/Data/DataRaw/METEO_VM_20XX_AL/"
plotdir = "/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/"
dataset = get_data(datadir)

dumptime = 60 #args.dumpt # in {10;30; 60} in min
 
# Raw data downloaded in format [unit/10min]

if idsensoreint in IdSensore[:,0]: #pluvio

    values = np.array(dataset['Valore Cumulato'])
    date = [x.to_pydatetime() for x in dataset['Data-Ora']]
    #print(date) ex. ...atetime.datetime(2023, 12, 30, 23, 40), datetime.datetime(2023, 12, 30, 23, 50), datetime.date...

    #idex_ok = np.where(values >= -1)[0]
    idex_nan = np.where(values <= -1)[0]
    #print('Warning : no data values \n', idex_nan)
    #values[idex_nan] = 0 #or use pandas fillna(0)
    #values = np.array(values)[idex_ok]
    #date = np.array(date)[idex_ok]

    dataset['Valore Cumulato'] = dataset['Valore Cumulato'].replace(to_replace=-999,value=0)
    values = np.array(dataset['Valore Cumulato'])
    
    fig, ax = plt.subplots(figsize=(10.0,7.0))
    ax.plot(date,values,label='RF')
    ax.plot(np.array(date)[idex_nan],np.array(values)[idex_nan],color='r',marker='.',linestyle='',label='Nodata')
    
    ax.set_ylim(ymin=0, ymax=20) #20 for 10min RF data
    #ax.set_xlabel('date')
    ax.set_ylabel('precipitation (mm)')
    ax.legend()
    ax.set_title('Cumulated precipitation in mm/10min - RawData - sensor '+idsensore)
    ax.grid(True)

    plotname = 'RF'+idsensore+'_10YA.png'
    fig.savefig(plotdir+plotname, bbox_inches='tight')

    # Correct NoData (minimal correction - error often last in long sequences)
    values[idex_nan]=(values[idex_nan-1]+values[idex_nan+1])/2

    if dumptime>10:
        dumpint = int(dumptime/10)
        nbval = int(len(date)/dumpint)
        dumpdate = [date[i] for i in range (dumpint,len(date),dumpint)]
        dumpvalues = [sum(values[dumpint*i:dumpint*(i+1)-1]) for i in range (0,nbval)]
        date = dumpdate
        values = dumpvalues

    fig, ax = plt.subplots(figsize=(10.0,7.0))
    ax.plot(date,values)
    ax.set_ylim(ymin=0, ymax=35) #args.valmaxplot) # 35 for 1h
    #ax.set_xlabel('date')
    ax.set_ylabel('precipitation (mm)')
    ax.set_title('Cumulated precipitation in mm/'+ str(dumptime) + 'min - sensor '+idsensore)
    ax.grid(True)

    #plt.show()

    #plotdir = ""
    plotname = 'RF'+idsensore+'_60YA.png' #YA = YearAll
    fig.savefig(plotdir+plotname, bbox_inches='tight')

    ### Annalysis by year/months

    #dataset['Valore Cumulato'] = dataset['Valore Cumulato'].replace(to_replace=-999,value=0)

    #dumpint = 6
    #nbval = int(len(date)/dumpint)
    #Dt = [date[i] for i in range (dumpint,len(date),dumpint)]
    #VC = [sum(values[dumpint*i:dumpint*(i+1)-1]) for i in range (0,nbval)]
    list_of_tuples = list(zip(date,values))
    weather = pd.DataFrame(list_of_tuples,columns=['Data-Ora','Valore Cumulato'])
    #dataset['Data-Ora'] = [date[i] for i in range (dumpint,len(date),dumpint)]
    #dataset['Valore Cumulato'] = [sum(values[dumpint*i:dumpint*(i+1)-1]) for i in range (0,nbval)]
    
    rainfall_by_year = split_data_by_year(weather)
    #print(rainfall_by_year)
    
    fig = plt.figure(figsize=(10.0, 7.0), dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    for year, annual_data in rainfall_by_year:
        days_since_jan1 = annual_data['days_since_jan1'].to_list()
        rainfall = annual_data['rainfall'].to_list()
        ax.plot(days_since_jan1, rainfall, label=year)
    
    # label the axes and title the plot
    ax.set_xlabel('days since jan 1st')
    ax.set_ylabel('precipitation (mm)')
    ax.set_ylim((-0.1,35.0))
    ax.set_title('Valmalenco /10min precipitation since Jan 1st')
    ax.grid(True)
    ax.legend()
    
    plotname = 'RF'+idsensore+'_60YC.png' #YearComparison
    fig.savefig(plotdir+plotname, bbox_inches='tight')
    #plt.show()

    ### Cumulative
    fig = plt.figure(figsize=(10.0, 7.0), dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    for year, annual_data in rainfall_by_year:
        days_since_jan1 = annual_data['days_since_jan1'].to_list()
        rainfall = annual_data['rainfall'].cumsum().to_list()
        ax.plot(days_since_jan1, rainfall, label=year)
    
    # label the axes and title the plot
    ax.set_xlabel('days since jan 1st')
    ax.set_xlim((-1,366))
    ax.set_ylabel('cumulative precipitation (mm)')
    ax.set_ylim((0,1200))
    ax.set_title('Valmalenco cumulative annual precipitation')
    ax.grid(True)
    ax.legend()
    
    plotname = 'RF'+idsensore+'_CYC.png'
    fig.savefig(plotdir+plotname, bbox_inches='tight')

    plt.show()

if idsensoreint in IdSensore[:,7]: #idro

    lab = 'WL'
    values = np.array(dataset[' Medio'])
    date = [x.to_pydatetime() for x in dataset['Data-Ora']]

    #idex_ok = np.where(values >= -1)[0]
    idex_nan = np.where(values <= -1)[0]
    #print('Warning : no data values \n', idex_nan)
    #values[idex_nan] = 0 #or use pandas fillna(0)
    #values = np.array(values)[idex_ok]
    #date = np.array(date)[idex_ok]

    dataset[' Medio'] = dataset[' Medio'].replace(to_replace=-999,value=np.nan)
    values = np.array(dataset[' Medio'])

    fig, ax = plt.subplots(figsize=(10.0,7.0))
    ax.plot(date,values,label=lab)
    ax.plot(np.array(date)[idex_nan],np.zeros(len(idex_nan)),color='r',marker='.',linestyle='',label='Nodata')
    
    ax.set_ylim(ymin=0, ymax=600)
    #ax.set_xlabel('date')
    ax.set_ylabel('water level (cm)')
    ax.legend()
    ax.set_title('Average water level in cm/10min - RawData - sensor '+idsensore)
    ax.grid(True)

    plotname = lab+idsensore+'_10YA.png'
    fig.savefig(plotdir+plotname, bbox_inches='tight')
    plt.show()

    # Correct NoData (minimal correction - error often last in long sequences)
    #values[idex_nan]=(values[idex_nan-1]+values[idex_nan+1])/2

    if (idsensoreint==16):

        # SIDRO
        # scala di deflusso - Sondrio Eifell
        # h from idsensore 16
        Qs = lambda h : 49.728*(h*1e-2-0.056)**1.67
        #Range misure di portata 0.247 < h(mm) < 0.785 :: in m ???
        #validita da 2013/09
        #ultimo aggiornamento 2019/09

        # rectangular SEZIONE with constant 2m/s velocity
        # 39-17
        Qr = lambda h : 22*h*1.e-2*1

        fig, ax = plt.subplots(figsize=(10.0,7.0))
        ax.plot(date,Qr(values),label='Qrect (v=1)')
        ax.plot(date,Qs(values),label='Qsidro')
        ax.plot(np.array(date)[idex_nan],np.zeros(len(idex_nan)),color='r',marker='.',linestyle='',label='Nodata')
        
        ax.set_ylim(ymin=0, ymax=150)
        #ax.set_xlabel('date')
        ax.set_ylabel('discharge (m3/s)')
        ax.legend()
        ax.set_title('Average discharge in m^3/s (/10min) - RawData - sensor '+idsensore)
        ax.grid(True)
        
        plotname = lab+idsensore+'_10YA_Q.png'
        fig.savefig(plotdir+plotname, bbox_inches='tight')
        plt.show()

