#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:09:58 2017

@author: pmarin
"""
import glob
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta

### Read data from the NPF Field Campaign
intN_pairs =  [(5,12),
              (7,30),
              (12,30),
              (30,140)]

filenames = glob.glob('/Volumes/Data/Act_Spectrum/NPF/*.xls')
aaaa_diam = []
for k in np.arange(0,len(filenames),1):
#for k in np.arange(0,4,1):
    
    fname = filenames[k] ; print(fname)
    
    df_temp = pd.read_excel(filenames[k],skiprows=1)   
    id_dN1 = df_temp.columns.get_loc("dN/dlogDp1")

    nbins = int((np.shape(df_temp)[1]-5)/2)
    id_dmin = df_temp.columns.get_loc("Dp1")
    id_dmax = df_temp.columns.get_loc("Dp"+str(nbins))
    id_dNmin = df_temp.columns.get_loc("dN/dlogDp1")
    id_dNmax = df_temp.columns.get_loc("dN/dlogDp"+str(nbins))
        
    if (id_dmax-id_dmin) != (id_dNmax-id_dNmin):
        print('Something Wrong in the DataFrame Ordering!')
        break

#### Test Plots        
#    fig = plt.figure(figsize=(10,7))
#    plt.contourf(df_temp['Data Time'].values,df_temp.iloc[0,2:2+59].values,np.transpose(df_temp.iloc[:,61:61+59].values),100)
#    #plt.contourf(df_temp['Data Time'].values,df_temp.iloc[0,2:2+59].values,np.log10(np.transpose(df_temp.iloc[:,61:61+59].values)),100)
#    plt.yscale('log')
#    plt.colorbar()
#    

#### Plot NPF Data ontop of Size distributions from Read_Temp_Variability_vFinal file        
#    factor = np.log(10)
#    for j in np.arange(0,np.shape(df_temp)[0]-5,24):       
#        date_string = fname[31:38]
#        cur_date = datetime.combine(datetime(int('20'+date_string[0:2]),int(date_string[2:4]),int(date_string[4:6])),df_temp.iloc[j,1]) + timedelta(hours=5)
#        
#        id_dma = np.abs(np.array(gtimes)-cur_date).argmin()
#        
##        fig = plt.figure(figsize=[10,6])
##        for i in np.arange(-4,5,1):   
##            plt.plot(df_temp.iloc[j+i,2:2+59].values,df_temp.iloc[j+i,61:61+59],'-b')
##        plt.plot(np.nanmean(df_temp.iloc[j-4:j+5,2:2+59].values,axis=0),np.nanmean(df_temp.iloc[j-4:j+5,61:61+59].values,axis=0),'-g',lw=4,label='NPF')
##        plt.plot(DP_full*1000,gdist_extrap_save[id_dma,:]*factor,c='0.40',lw=7,label='Step 1: Extrapolation')
##        plt.plot(DP_full*1000,gdist_extrap_cpc_save[id_dma,:]*factor,c='0.55',lw=5,label='Step 2: CPC Efficiency Scaling')
##        plt.plot(DP_full*1000,gdist_extrap_cpc_last_save[id_dma,:]*factor,c='0.70',lw=3,label='Step 3: Remaining Difference: '+str(round(cur_cpcint_dif[id_dma]))+' / '+str(round(cur_cpcint_pdif[id_dma]*100))+'% Particles')
##        plt.plot(DP_full*1000,gdist_extrap_cpc_last_est_save[id_dma,:]*factor,'--',c='r',lw=2,label='Step 4: CPC Efficiency Scaling')
##        plt.plot(standDP*1000,gdist[id_dma,:]*factor,'-k',lw=2,label='DMA_unadjusted')
##        #plt.plot(DP_full*1000,gdist_extrap_cpc_est_smooth[id_dma,:],c='r',lw=1,label='Step 5: Smooth')
##
##        #plt.contourf(df_temp['Data Time'].values,df_temp.iloc[0,2:2+59].values,np.log10(np.transpose(df_temp.iloc[:,61:61+59].values)),100)
##        plt.legend(fontsize=10)
##        plt.xscale('log')
##        plt.title(str(cur_date)+' | '+str(round(cur_cpcint_dif[id_dma])))
##        plt.savefig('/Users/pmarin/Documents/ACT_Spectrum_Test/v0605_Size_Dist_'+str(cur_date)+'_EZRA_addpart_expon4_2W_curvefitDP2_windividuals_MeanCPCDMA_withCCN.png')
##        plt.close(fig)
##        
    intN_npf = np.zeros((np.shape(df_temp)[0],np.shape(intN_pairs)[0]))
    intN_all = np.zeros(np.shape(df_temp)[0])
    cur_date = []    
    for j in np.arange(0,np.shape(df_temp)[0]):       

        #date_string = str(df_temp['filename'][j])
        date_string = fname[31:38]
        cur_date.append(datetime.combine(datetime(int('20'+date_string[0:2]),int(date_string[2:4]),int(date_string[4:6])),df_temp.iloc[j,1]) + timedelta(hours=5))
        diams = np.asarray(df_temp.iloc[j,id_dmin:id_dmax].values).astype(np.float)
        diams_edge = (diams[1:]+diams[:-1])/2
        diams_edge = np.hstack((diams_edge[0]-(diams_edge[1]-diams_edge[0]),diams_edge,diams_edge[-1]+(diams_edge[-1]-diams_edge[-2])))
        dlog_diams = np.diff(np.log10(diams_edge))

        for i in np.arange(0,np.shape(intN_pairs)[0]):
            id_npf1 = np.abs(diams-intN_pairs[i][0]).argmin()
            id_npf2 = np.abs(diams-intN_pairs[i][1]).argmin()
            #print(diams[id_npf1],diams[id_npf2])
 
            if df_temp['total concentration (#/cm3)'][j] > 1000000:
                intN_npf[j,i] = np.nan
            else:
                intN_npf[j,i] = np.nansum(df_temp.iloc[j,id_dN1+id_npf1:id_dN1+id_npf2].values*dlog_diams[id_npf1:id_npf2])

        intN_all[j] = np.nansum(df_temp.iloc[j,id_dN1:id_dN1+len(diams)].values*dlog_diams)

######## PLOT NPF Size Distribution for each date        
#        fig = plt.figure(figsize=(10,7))
#        plt.plot(diams,df_temp.iloc[j,id_dN1:id_dN1+len(diams)],'-k',lw=4)
#        plt.xscale('log')
#        plt.xlabel('Diameter (micron)')
#        plt.ylabel('dN/dlog(Dp)')
#        plt.savefig('/Volumes/Data/Act_Spectrum/NPF/Plots/NPF_Size_Distribution_'+str(cur_date[j])+'.png')
#        plt.close(fig)
        
    # add new columns to dataframe with integrated number and datetime variables
    df_temp['time_adj'] = cur_date
    for i in np.arange(0,np.shape(intN_pairs)[0]):
        df_temp['nm'+str(intN_pairs[i][0])+str(intN_pairs[i][1])] = intN_npf[:,i]
    df_temp['all'] = intN_all

    if k == 0:
        df_npf = df_temp;
    else:
        df_npf = pd.concat([df_npf,df_temp],axis=0)
        
    del(df_temp)

df_npf = df_npf.set_index('time_adj')    

date_npf = []
start_date = datetime(2013,4,19,0,0,0)
for i in np.arange(0,len(filenames)*288,1):
    dt = np.int(i*5)
    date_npf.append(start_date+timedelta(minutes=dt))

# Create dataframe that is used in Read_Temp_Variability_vFinal       
dfn_120M = df_npf.resample('120Min',loffset='60Min').mean()
new_date_index_2H = [datetime(2013,4,19,1,0,0) + timedelta(hours=2*x) for x in range(0,500)]
dfn_120M = dfn_120M.reindex(new_date_index_2H, fill_value = np.nan)
