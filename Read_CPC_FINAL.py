#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:09:58 2017

@author: pmarin
"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from datetime import timedelta
import pandas as pd
import itertools
from scipy.io import netcdf
import netCDF4 as nc
from netCDF4 import Dataset
import copy

######################################################################
# CPC DATA - sgpnoaaaosC1.b1 datastream
######################################################################

pathname3 = '/Volumes/Data/Act_Spectrum/CPC/'  # May 2007 - Dec. 2013
#pathname3 = '/Volumes/Data/Act_Spectrum/CPC_Other/'
#g = MFDataset(pathname2+"sgpmet*cdf")
g = os.listdir(pathname3)
c_time1 = list(); ccn1 = list(); cpc1 = list(); flags1 = list()
for i in g:
    print('Reading Filename:',i)
    if i.startswith('sgp'):
    #if i.startswith('sgpnoaaaosC1.b1.20090505'):
    #if i.startswith('sgpnoaaaosC1.b1.20121125'):
              
        fi = Dataset(pathname3+i, "r", format="NETCDF4")       
        #nc2 = netcdf.netcdf_file(pathname3+i,'r')

        #asfdgsd

        base = nc.num2date(fi.variables['base_time'][:], units=fi.variables['base_time'].units)
        to_temp = fi.variables['time_offset'][:]    
        arr = np.array([base + timedelta(seconds=i) for i in to_temp])
        for item in arr:
            c_time1.append(item)    
        ccn1.append(fi.variables['N_CCN_1'][:])
        cpc1.append(fi.variables['N_CPC_1'][:])
        #flags1.append(fi.variables['flags_CMDL'][:])

        fi.close(); #nc2.close()

ccn1 = np.asarray(list(itertools.chain(*ccn1))); ccn1[ccn1<=0] = np.nan
cpc1 = np.asarray(list(itertools.chain(*cpc1))); cpc1[cpc1<=0] = np.nan
cpc1[cpc1<=200] = np.nan
cpc1[cpc1>=75000] = np.nan

# Fix the switch in datastreems during this time period
id_fix1 = np.abs(np.array(c_time1) - datetime(2009,4,28,19,30)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2009,5,15,18,31)).argmin()
for i in np.arange(id_fix1,id_fix2):
    temp1 = ccn1[i]; temp2 = cpc1[i];
    ccn1[i] = temp2; cpc1[i]=temp1

cpc1_fix = copy.deepcopy(cpc1)

# Fix the switch in leakage in the CPC from -- DQR Report D120608.1
id_fix1 = np.abs(np.array(c_time1) - datetime(2010,10,1,15,40)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2011,4,18,20,10)).argmin()
for i in np.arange(id_fix1,id_fix2):
    #cpc1_fix[i] = cpc1_fix[i] / 0.8;
    cpc1_fix[i] = np.nan;

# Fix clearly low values in the CPC during this time period 
# (both sharp deceases and rises at the beginning and end of this period)
id_fix1 = np.abs(np.array(c_time1) - datetime(2012,11,24,13,0)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2013,1,31,17,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc1_fix[i] = np.nan;

# Fix clearly low values in the CPC during this time period 
# (both sharp deceases and rises at the beginning and end of this period)
# CPC values less than the CCN data
id_fix1 = np.abs(np.array(c_time1) - datetime(2010,4,4,0,0)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2010,4,6,21,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc1_fix[i] = np.nan;

# Fix clearly low values in the CPC during this time period 
# (both sharp deceases and rises at the beginning and end of this period)
# CPC values less than the CCN data
id_fix1 = np.abs(np.array(c_time1) - datetime(2010,4,10,0,0)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2010,4,27,20,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc1_fix[i] = np.nan;

# Fix clearly low values in the CPC during this time period 
# (both sharp deceases and rises at the beginning and end of this period)
# CPC values less than the CCN data
id_fix1 = np.abs(np.array(c_time1) - datetime(2010,6,30,23,0)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2010,7,14,15,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc1_fix[i] = np.nan;

# Fix clearly low values in the CPC during this time period 
# (both sharp deceases and rises at the beginning and end of this period)
# CPC values less than the CCN data
id_fix1 = np.abs(np.array(c_time1) - datetime(2010,9,27,21,0)).argmin()
id_fix2 = np.abs(np.array(c_time1) - datetime(2010,9,28,16,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc1_fix[i] = np.nan;

#
#id_fix1 = np.abs(np.array(c_time1) - datetime(2012,11,24,13,0)).argmin()
#id_fix2 = np.abs(np.array(c_time1) - datetime(2013,1,31,17,0)).argmin()
#fig = plt.figure(figsize=(20,5));
#plt.plot(c_time1[id_fix1-100:id_fix2+100],cpc1[id_fix1-100:id_fix2+100],label='CPC')
#plt.plot(c_time1[id_fix1-100:id_fix2+100],ccn1[id_fix1-100:id_fix2+100],label='CCN')
#plt.plot(c_time1[id_fix1-100:id_fix2+100],cpc1_fix[id_fix1-100:id_fix2+100],label='CPC_fix')
##plt.plot(c_time2,cpc2,label='CPC2')
##plt.plot(c_time3,cpc3,label='CPC3')
##plt.plot(times[1:25650],intN[1:25650],label='DMA')
#plt.legend()
#plt.xlim([datetime(2010,4,1),datetime(2010,9,30)])
#plt.ylim([0,10000])
##
##
#fig = plt.figure(figsize=(10,5));
#plt.plot(c_time1,cpc1_fix/ccn1)
##plt.plot(c_time1,cpc1/ccn1)
#plt.plot(c_time1,np.zeros(len(cpc1))+1,'-k')
#plt.ylim([0.1,100])
#plt.yscale('log')
#plt.xlim([datetime(2013,1,1),datetime(2014,1,1)])
#

#fig,ax = plt.subplots(4,1,figsize=(18,12));
#cnt = 0
#for i in np.arange(2009,2012):
#    ax[cnt].plot(c_time1,cpc1,label='CPC1')
#    ax[cnt].plot(c_time1,ccn1,label='CCN1')
#    #ax[0].plot(c_time2,cpc2,label='CPC2')
#    ax[cnt].plot(c_time1,cpc1_smooth[29:len(cpc1_smooth)-30],label='CPC1_sm')
#    ax[cnt].plot(np.array(times[1:]),intN[1:],label='DMA')
#    ax[cnt].plot(np.array(times[1:]),ccn_max[1:],label='DMA-CCN')
#    ax[cnt].plot(np.array(times[1:]),cpc_data[1:,0],label='DMA-CPC')
#    ax[cnt].legend()
#    ax[cnt].set_xlim([datetime(i,10,1),datetime(i+1,4,1)])
#    ax[cnt].set_ylim([0,20000])
#    cnt = cnt + 1



######################################################################
# CPC DATA - sgpaosC1.b1 datastream 
######################################################################
pathname3 = '/Volumes/Data/Act_Spectrum/CPC_C1_b1/'   # March 2011 - Dec. 2017

g = os.listdir(pathname3)
c_time2 = list(); cpc2 = list(); flags2 = list();  f2_flags2 = list();
for i in g:
    print('Reading Filename:',i)
    if i.startswith('sgp'):
              
        fi = Dataset(pathname3+i, "r", format="NETCDF4")       
        #nc2 = netcdf.netcdf_file(pathname3+i,'r')

        #dsfgsgfdsg

        base = nc.num2date(fi.variables['base_time'][:], units=fi.variables['base_time'].units)
        to_temp = fi.variables['time_offset'][:]    
        arr = np.array([base + timedelta(seconds=i) for i in to_temp])
        for item in arr:
            c_time2.append(item)    
        cpc2.append(fi.variables['concentration'][:])
        flags2.append(fi.variables['qc_concentration'][:])
        f2_flags2.append(fi.variables['f2_flags'][:])
        #dcf2.append(nc2.variables['dilution_correction_factor'].data)
        #cpc2_flow.append(nc2.variables['cpc_flow'].data)

        fi.close(); #nc2.close()

cpc2 = np.asarray(list(itertools.chain(*cpc2))); 
cpc2[cpc2<=0] = np.nan
cpc2[cpc2>=75000] = np.nan
cpc2[cpc2<=200] = np.nan

cpc2_fix = copy.deepcopy(cpc2)

# Fix clearly low values in the CPC during this time period 
# (both sharp deceases and rises at the beginning and end of this period)
id_fix1 = np.abs(np.array(c_time2) - datetime(2012,11,24,13,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2013,1,31,17,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2013,7,25,0,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2013,7,25,4,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2013,7,25,13,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2013,7,29,18,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2013,10,4,15,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2013,10,4,18,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2013,7,25,12,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2013,7,29,18,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2013,8,9,9,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2013,8,12,18,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2014,6,4,2,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2014,6,5,22,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2014,6,16,12,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2014,6,17,15,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2014,6,21,1,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2014,7,17,20,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

id_fix1 = np.abs(np.array(c_time2) - datetime(2015,6,22,0,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2015,6,29,12,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

# Error report (D160619.7 --- particle counts lower than other instruments -- probably due to a combination of contaminated butanol and a dirty focusing nozzle)
id_fix1 = np.abs(np.array(c_time2) - datetime(2016,5,8,12,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2016,5,17,21,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

# Error report (D160707.4 --- particle counts lower due to a weak laser)
id_fix1 = np.abs(np.array(c_time2) - datetime(2016,6,30,20,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2016,7,18,16,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

# Error report (D170619.11 --- particle counts lower than other instruments -- probably due to a combination of contaminated butanol and a dirty focusing nozzle)
id_fix1 = np.abs(np.array(c_time2) - datetime(2017,6,16,0,0)).argmin()
id_fix2 = np.abs(np.array(c_time2) - datetime(2017,6,16,16,0)).argmin()
for i in np.arange(id_fix1,id_fix2):
    cpc2_fix[i] = np.nan;

flags2 = np.asarray(list(itertools.chain(*flags2)))
flag_bin = list()
for i in np.arange(0,len(flags2)):
    flag_bin.append(format(flags2[i],'#010b'))

flag_bin2 = np.zeros((8,len(flags2)))
for i in np.arange(0,len(flags2)):
    flag_bin2[7,i] = int(flag_bin[i][2:3])
    flag_bin2[6,i] = int(flag_bin[i][3:4])
    flag_bin2[5,i] = int(flag_bin[i][4:5])
    flag_bin2[4,i] = int(flag_bin[i][5:6])
    flag_bin2[3,i] = int(flag_bin[i][6:7])
    flag_bin2[2,i] = int(flag_bin[i][7:8])
    flag_bin2[1,i] = int(flag_bin[i][8:9])
    flag_bin2[0,i] = int(flag_bin[i][9:10])

# Eliminate data due to vacuum error, instrument not ready or missing values
for i in np.arange(0,len(cpc2)):
    if (flag_bin2[7,i] == 1) or (flag_bin2[6,i] == 1) or (flag_bin2[1,i] == 1):
        cpc2_fix[i] = np.nan

#                  min_concentration_warning: 100.0
#                  max_concentration_warning: 8000.0
#                  negative_delta_warning: -600.0
#                  flag_method: bit
#                  bit_1_description: Value is equal to missing_value.
#                  bit_1_assessment: Bad
#                  bit_2_description: Value is less than the valid_min.
#                  bit_2_assessment: Bad
#                  bit_3_description: Value is greater than the valid_max.
#                  bit_3_assessment: Bad
#                  bit_4_description: Difference between current and previous sample <= negative_delta_warning
#                  bit_4_assessment: Indeterminate
#                  bit_4_comment: Check not performed if negative_delta_warning is -9999
#                  bit_5_description: Concentration is less than min_concentration_warning
#                  bit_5_assessment: Indeterminate
#                  bit_6_description: Concentration is greater than max_concentration_warning
#                  bit_6_assessment: Indeterminate
#                  bit_7_description: f2_flags indicates "not_ready"
#                  bit_7_assessment: Bad
#                  bit_8_description: f2_flags field indicates "vacuum_error"
#                  bit_8_assessment: Bad



#fig = plt.figure(figsize=(20,5));
##plt.plot(c_time1,cpc1,label='CPC')
##plt.plot(c_time1,ccn1,label='CCN')
##plt.plot(c_time1,cpc1_fix,label='CPC_fix')
##plt.plot(c_time2,ccn2,label='CPC2')
#plt.plot(c_time2,cpc2,label='CPC2')
#plt.plot(c_time2,cpc2_fix,label='CPC2_fix')
##plt.plot(c_time3,cpc3,label='CPC3')
##plt.plot(times[1:25650],intN[1:25650],label='DMA')
#plt.legend()
##plt.xlim([datetime(2010,4,1),datetime(2010,9,30)])
#plt.ylim([0,10000])

#f2_flags2 = np.asarray(list(itertools.chain(*f2_flags2)))
#f2_flag_bin = list()
#for i in np.arange(0,len(f2_flags2)):
#    f2_flag_bin.append(format(f2_flags2[i],'#010b'))
#
#f2_flag_bin2 = np.zeros((8,len(f2_flags2)))
#for i in np.arange(0,len(flags2)):
#    f2_flag_bin2[7,i] = int(f2_flag_bin[i][2:3])
#    f2_flag_bin2[6,i] = int(f2_flag_bin[i][3:4])
#    f2_flag_bin2[5,i] = int(f2_flag_bin[i][4:5])
#    f2_flag_bin2[4,i] = int(f2_flag_bin[i][5:6])
#    f2_flag_bin2[3,i] = int(f2_flag_bin[i][6:7])
#    f2_flag_bin2[2,i] = int(f2_flag_bin[i][7:8])
#    f2_flag_bin2[1,i] = int(f2_flag_bin[i][8:9])
#    f2_flag_bin2[0,i] = int(f2_flag_bin[i][9:10])

#id1 = np.abs(np.array(c_time2) - datetime(2015,7,1)).argmin()
#id2 = np.abs(np.array(c_time2) - datetime(2016,1,1)).argmin()
##
#fig,ax = plt.subplots(2,1,sharex=True);
#for i in np.arange(0,8):
#    ax[0].plot(c_time2[id1:id2],flag_bin2[i,id1:id2]*(i+1),'o',label=str(i+1))
#ax[0].legend()
#ax[1].plot(c_time2[id1:id2],cpc2[id1:id2],label='CPC2')
#ax[1].plot(c_time2[id1:id2],cpc2_fix[id1:id2],label='CPC2_FIX')
#ax[1].legend()
#ax[1].set_yscale('log')
#ax[1].set_ylim([1,100000])

######################################################################
# CPC DATA - sgpnoaaaosavgC1.b1 datastream (avergaged / smoothed fields)
######################################################################

#pathname3 = '/Volumes/Data/AOS_AVG3/'
##pathname3 = '/Volumes/Data/Act_Spectrum/CPC_Other/'
##g = MFDataset(pathname2+"sgpmet*cdf")
#g = os.listdir(pathname3)
#c_time3 = list(); ccn3 = list(); cpc3 = list(); flags3 = list()
#for i in g:
#    if i.startswith('sgp'):
#        if i.endswith('cdf') or i.endswith('nc'):
#        #if i.startswith('sgpnoaaaosC1.b1.20090505'):
#        #if i.startswith('sgpnoaaaosC1.b1.20121125'):
#            print('Reading Filename:',i)              
#            fi = Dataset(pathname3+i, "r", format="NETCDF4")       
#            nc2 = netcdf.netcdf_file(pathname3+i,'r')
#            base = nc.num2date(fi.variables['base_time'][:], units=fi.variables['base_time'].units)
#            to_temp = fi.variables['time_offset'][:]    
#            arr = np.array([base + timedelta(seconds=i) for i in to_temp])
#            for item in arr:
#                c_time3.append(item)    
#            #ccn3.append(nc2.variables['N_CCN_1'].data)
#            cpc3.append(nc2.variables['N_CPC_1'].data)
#            #flags3.append(nc2.variables['flags_CMDL'][:])
#
#            fi.close(); nc2.close()
#
##ccn3 = np.asarray(list(itertools.chain(*ccn3))); ccn3[ccn3<=0] = np.nan
#cpc3 = np.asarray(list(itertools.chain(*cpc3))); cpc3[cpc3<=0] = np.nan
#
##fig = plt.figure(figsize=(10,5));
##plt.plot(c_time3,cpc3,label='CPC3')
#plt.legend()
#plt.xlim([datetime(2012,10,1),datetime(2013,3,30)])
#plt.ylim([0,10000])


######################################################################
# Combine Datastreams
######################################################################
# Convert data to pandas dataframe
df1_cpc = pd.DataFrame.from_dict({'time': c_time1,
                                 'cpc': cpc1_fix,
                                 'ccn': ccn1})                                 
df1_cpc = df1_cpc.set_index('time')

df2_cpc = pd.DataFrame.from_dict({'time': c_time2,
                                 'cpc2': cpc2_fix})                                 
df2_cpc = df2_cpc.set_index('time')

df_cpc = pd.concat([df1_cpc,df2_cpc],axis=0)
df_cpc = df_cpc.resample(rule='5Min',how=['mean'])
df_cpc['cpc_avg'] = np.nanmax([df_cpc.cpc.values,df_cpc.cpc2.values],axis=0)

dfc_1440M = df_cpc.resample('1440Min',loffset='720Min').mean()
from datetime import datetime, timedelta
new_date_index_24H = [datetime(2007,1,1,12,0,0) + timedelta(hours=24*x) for x in range(0,3830)]
dfc_1440M = dfc_1440M.reindex(new_date_index_24H, fill_value = np.nan)
dfc_1440M['yr'] = dfc_1440M.index.year; dfc_1440M['mo'] = dfc_1440M.index.month; dfc_1440M['hr'] = dfc_1440M.index.hour

import pickle
pickle.dump(dfc_1440M,open('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/df_1440M_CPC_FINAL.p',"wb"))

#
fig = plt.figure(figsize=(20,5));
plt.plot(df_cpc.index,df_cpc.cpc_avg,label='CPC_AVG')
plt.plot(df_cpc.index,df_cpc.cpc,'o',label='CPC1')
plt.plot(df_cpc.index,df_cpc.cpc2,'x',label='CPC2')
plt.legend()

fig = plt.figure(figsize=(20,5));
plt.plot(dfc_1440M.index,dfc_1440M.cpc_avg,label='CPC_AVG')
plt.plot(dfc_1440M.index,dfc_1440M.cpc,'o',label='CPC1')
plt.plot(dfc_1440M.index,dfc_1440M.cpc2,'x',label='CPC2')
plt.legend()