#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 10:46:02 2018

@author: pmarin
"""

import datetime

plt.rcParams.update(plt.rcParamsDefault)

id_period = 1
pr_desired = 0.5 #days

# First Frequency of Interst
find_freq = np.abs(pr_dys-pr_desired).argmin() + 1
print(str(pr_dys[find_freq-1])+' day cycle')
hours_for_analysis = pr_dys[find_freq-1]*24;

# Second frequency of interest
pr_desired_2 = 1.0 #days
find_freq_2 = np.abs(pr_dys-pr_desired_2).argmin() + 1
print(str(pr_dys[find_freq_2-1])+' day cycle')

# Load data from other function output
PWR = pwr_s_comb[id_period]
IFFT_RE = ifft_re_s_comb[id_period]
TIMES = chunktime_s_comb[id_period]

# Find powerful times in association with the first frequency of interest
pwr_thresh = copy.deepcopy(F99red2[id_period][find_freq-1])

pwr_thresh = 0.10
# Find times when there is high power
cnt_pwr = 0
id_pwr = []
for i in np.arange(0,np.shape(TIMES)[0]):   
    if (PWR[i,find_freq_2-1] > pwr_thresh):
        cnt_pwr = cnt_pwr+1
        id_pwr.append(i)
id_pwr = np.asarray(id_pwr)

# Calculate sum of all the different frequencies (Make them equal to periods)
IFFT_SUM2 = 2*np.nansum(IFFT_RE,axis=0)
for i in np.arange(0,np.shape(IFFT_SUM2)[0]):
    for j in np.arange(0,np.shape(IFFT_SUM2)[1]):
        IFFT_SUM2[i,j] = IFFT_SUM2[i,j] - IFFT_RE[0,i,j] - IFFT_RE[42,i,j]

# Calculate sum of all the different frequencies (Make them equal to periods)
fig,ax = plt.subplots(1,1,figsize=(20,6))
plt.plot(pwr_s_comb[id_period][id_pwr,find_freq-1])
plt.plot([0,cnt_pwr],[pwr_thresh,pwr_thresh])
  
ids_time = np.arange(0,84)
            
import matplotlib.dates as mdates
#for i in np.arange(0,len(id_pwr),1):
for i in np.arange(0,len(id_pwr)):
    
    print(pwr_s_comb[id_period][id_pwr[i],find_freq-1])
    datestr = str(TIMES[id_pwr[i]][0])[0:10]
    YY = str(TIMES[id_pwr[i]][0])[0:4]
    MM = str(TIMES[id_pwr[i]][0])[5:7]
    print(datestr)
    

    fig,ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(pr_dys*24,PWR[id_pwr[i],:])
    ax.plot(pr_dys*24,F99red2[id_period])
    ax.set_xscale('log')
    ax.set_xlabel('Period')
    ax.set_ylabel('Norm. Power')

    #ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax.set_xticks([4,6,12,24,48,84])
    ax.set_xticklabels([4,6,12,24,48,84])   
    ax.set_ylim([0.001,0.2])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim([dt*2/24*days_conv * 0.97,num_days*days_conv/2 * 1.03])
    plt.tight_layout()
#    plt.savefig('/Users/pmarin/Aerosol_Power_Spectra_Breakdown/Power_Spectra_'+datestr+'.png')
#    plt.close(fig)
#    
    
    for num_freq in np.arange(0,40):
        plt.rcParams.update({'font.size': 20})
        fig,ax = plt.subplots(1,1,figsize=(16,7))
    
        ax.plot(TIMES[id_pwr[i]][ids_time],IFFT_SUM2[id_pwr[i],ids_time]-IFFT_RE[0,id_pwr[i],:],'-',c='0.5',lw=10,label='Original Data')

        if num_freq > 0:
   
            now_PWR = PWR[id_pwr[i],:]
            now_PWR[0] = 0
            freq_largest_power = sorted(range(len(now_PWR)), key=lambda i: now_PWR[i])[-num_freq:]
            freq_largest_power = freq_largest_power[::-1]
    
            pwr_sum_10 = np.zeros(84)        
            for j in freq_largest_power:
                if j == 0:
                    temp_pwr_sum_10 = IFFT_RE[j+1,id_pwr[i],:]
                else:
                    temp_pwr_sum_10 = 2*IFFT_RE[j+1,id_pwr[i],:]            
                pwr_sum_10 = pwr_sum_10 + temp_pwr_sum_10 
            ax.plot(TIMES[id_pwr[i]][ids_time],pwr_sum_10[ids_time],'-',c='0.1',lw=8,label='Sum of Cycles')
    
            pwr_sum_10 = np.zeros(84)        
            for j in freq_largest_power:
                if j == 0:
                    temp_pwr_sum_10 = IFFT_RE[j+1,id_pwr[i],:]
                else:
                    temp_pwr_sum_10 = 2*IFFT_RE[j+1,id_pwr[i],:]
    
                ax.plot(TIMES[id_pwr[i]][ids_time],temp_pwr_sum_10[ids_time],'-',lw=2,label=str(round(per_hrs[j],1))+' Hr Cycle')#:'+str(round(PWR[id_pwr[i],j],3)))
                pwr_sum_10 = pwr_sum_10 + temp_pwr_sum_10 
# 
      
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
        ax.set_xlabel('Day | Time ('+MM+'-'+YY+')')
        ax.set_ylabel('Total Aerosol Concentration (N$_{T}$) \n Anomaly (# cm$^{-3}$)')
        ax.set_ylim([-7500,12500])
        ax.legend(bbox_to_anchor=(1.04,1))

        plt.tight_layout()
        plt.subplots_adjust(right=0.63)
#        plt.savefig('/Users/pmarin/Aerosol_Power_Spectra_Breakdown/'+datestr+'_Frequency_'+str(num_freq)+'.png')
#        plt.close(fig)
    

#    plt.savefig('/Volumes/Data/Act_Spectrum/Test_Sums/'+'Test_'+str(TIMES[id_pwr[i]][0])[0:10]+'_'+str(int(id_pwr[i]))+'.png')
#    plt.close(fig)
#
#
#
#    fig,ax = plt.subplots(1,1,figsize=(13,4.5))
#    
#    #plt.plot(TIMES[id_pwr[i]],IFFT_SUM[id_pwr[i],:]/np.hanning(84),'-k',lw=8,label='Original Data')
#    #ax.plot(TIMES[id_pwr[i]],IFFT_SUM2[id_pwr[i],:]/np.hanning(84),'-k',lw=5,label='Original Data')
#    ax.plot(TIMES[id_pwr[i]],IFFT_SUM2[id_pwr[i],:],'-',c='0.2',lw=7,label='Original Data + Hanning')
#    
##    temp_sum = IFFT_RE[0,id_pwr[i],:]
##    for j in np.arange(1,find_freq+1):
##        temp_sum = temp_sum + 2*IFFT_RE[j,id_pwr[i],:]
##    ax[0].plot(TIMES[id_pwr[i]],temp_sum,'-c',lw=6,label='0-14 Frequencies')
#    
#    ax.plot(TIMES[id_pwr[i]],2*IFFT_RE[find_freq,id_pwr[i],:]+2*IFFT_RE[find_freq_2,id_pwr[i],:],'-g',lw=6,label='12-hour + 24-hour cycles')
#    ax.plot(TIMES[id_pwr[i]],2*IFFT_RE[find_freq_2,id_pwr[i],:],'-y',lw=4,label='24-hour cycle')#: '+str(round(PWR[id_pwr[i],find_freq_2-1],3)))
#    ax.plot(TIMES[id_pwr[i]],2*IFFT_RE[find_freq,id_pwr[i],:],'-c',lw=4,label='12-hour cycle')#: '+str(round(PWR[id_pwr[i],find_freq-1],3)))
#    ax.legend()
#    ax.xaxis.set_ticks([datetime.datetime(2012,2,22,12) + datetime.timedelta(hours=int(i*12)) for i in np.arange(0,15)])
#
#    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
#
#    ax.set_xlabel('UTC Time (Feb. 2012 Day|Hour)')
#    ax.set_ylabel('N$_{7-30nm}$ Concentration Anamoly \n (# $cm^{-3}$)')
##    ax.set_title(str(TIMES[id_pwr[i]][0])[0:10])
#    ax.grid()
#    plt.tight_layout()
#    plt.savefig('/Volumes/Data/Act_Spectrum/PUB_Example_12+24hour_Cycle.png')
#    plt.savefig('/Volumes/Data/Act_Spectrum/PUB_Example_12+24hour_Cycle.eps')
#
#
#
#
#IFFT_SUM2 = 2*np.nansum(IFFT_RE,axis=0)
#for i in np.arange(0,np.shape(IFFT_SUM2)[0]):
#    for j in np.arange(0,np.shape(IFFT_SUM2)[1]):
#        IFFT_SUM2[i,j] = IFFT_SUM2[i,j] - IFFT_RE[0,i,j] - IFFT_RE[42,i,j]
#
#fig = plt.figure(figsize=(20,6))
#plt.plot(df_120M_2011.index,df_120M_2011.intN0.values-np.nanmean(df_120M_2011.intN0.values),'-k')
#for i in np.arange(0,len(id_pwr),1):
#    plt.plot(TIMES[id_pwr[i]],IFFT_SUM[id_pwr[i],:],'--r')
#    plt.plot(TIMES[id_pwr[i]],IFFT_SUM[id_pwr[i],:]/np.hanning(84),'--g')
#
#
#    
#       
#j_id = 11
#fig = plt.figure(figsize=(20,6))
#for i in np.arange(0,43):
#    plt.plot(TIMES[j_id],IFFT_RE[i,j_id,:])