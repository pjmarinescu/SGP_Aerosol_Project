#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:19:10 2018
@author: pmarin
"""

# Calculate Power Spectrum of Data
import math
import numpy as np
import scipy
import copy
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
  
def power_spectrum(df_in,dt,chunk_size,overlap):        

# FOR TESTING
##    df_in = dfa_120M_JJA_2011.total
#    df_in = df_120M_2012.intN_v2
#    dt = 2 # hours 
#    chunk_size = 7 # days
#    overlap = 0.50

#    df_in = cur_data
#    dt = 24
#    chunk_size = 28
#    overlap = 0.5


    N_period = (chunk_size*24/dt); # days * hours / dt  # lowest frequency, largest wavenumber
    if N_period % 2 == 0:
        T_period = N_period + 1 # T needs to be odd / N needs to be even
    else:
        T_period = N_period
        N_period = T_period - 1
    sizeChunks = int(N_period) # Number of data points to go into each chunk (limits the ability to assess low frequency signals)

    # Reformat df_in
    testdata = df_in.values
    testdata_time = df_in.index
    testdata_hr = df_in.index.hour.values
    testdata = (testdata-np.nanmean(testdata)) 

    lag1_pd = df_in.autocorr(1)

    # Do power spectrum on anomalies of data, assumption in 655 notes
    # also Hanning Window does wierd things if the mean is ~0
    ######################################################        
    # Input and calculate parameters for analysis
    ######################################################
    T_all = len(testdata) # data record length
    Hanning = np.hanning(sizeChunks); # hanning window for the size of desired chunks
    wave_num = np.arange(1,sizeChunks/2+1) 
    freq = wave_num / sizeChunks; # k / T, recall 1:N/2 are the wavenumbers
    per = 1/freq; per_hrs = 1/(freq/dt); per_days = per_hrs/24
    wave_num_Nyq = (sizeChunks-1)/2 # 
    freq_Nyq = 1/(2*dt)  # Nyquist Frequency #/hr OR samplying frequency of data / 2

    # Determine how many chunks for a given overlap %
    nChunks = math.floor(T_all/sizeChunks/(1-overlap))
    idx = np.zeros((nChunks+1,2)); 
    for n in np.arange(0,nChunks+1):
        # Get data from each chunk
        if n == 0:
            idx[n,:] = [0, sizeChunks] # First chunk will always be 0 to sizeChunk            
        else:
            idx[n,:] = [int(idx[n-1,0]+((1-overlap)*sizeChunks)), int(idx[n-1,0]+((1-overlap)*sizeChunks)) + sizeChunks]   # No overlap
        
        if idx[n,1] > T_all:
            break
    nChunks = n

    #####################################################
    # Specifiy variables needed to each dataset for chunking
    #####################################################   
    cnt = 0;      
    idx = np.zeros((nChunks,2)); 
    powers = np.zeros((nChunks,sizeChunks));
    power_mod = np.zeros((nChunks,sizeChunks)); 
    fft_data = np.zeros((nChunks,sizeChunks))
    fft_data_time = np.zeros((nChunks,sizeChunks))
    fft_data_hr = np.zeros((nChunks,sizeChunks))
    power_var = np.zeros((nChunks,int(sizeChunks/2)))
    power_var_app = np.zeros((1,int(sizeChunks/2)))
    # For each chunk of data, save the 
    ifft_data_re = np.zeros((int(sizeChunks/2+1),nChunks,sizeChunks)) # num frequencies, num chunks, size of each chunk
    ifft_data_im = np.zeros((int(sizeChunks/2+1),nChunks,sizeChunks)) # num frequencies, num chunks, size of each chunk
    pk_time = np.zeros((2,int(sizeChunks/2+1),nChunks))

    Te = np.zeros(nChunks); lagTe = np.zeros(nChunks); lag1 = np.zeros(nChunks);

    chunk_datetime = []
    # Separate each dataset into its Chunks, and calculate power spectrum
    for n in np.arange(0,nChunks):
        # Get data from each chunk
        if n == 0:
            idx[n,:] = [0, sizeChunks]   # No overlap            
        else:
            idx[n,:] = [int(idx[n-1,0]+((1-overlap)*sizeChunks)), int(idx[n-1,0]+((1-overlap)*sizeChunks)) + sizeChunks]   # No overlap

        testdata_chunk = copy.deepcopy(testdata[int(idx[n][0]):int(idx[n][1])])
        fft_data[n,:] = testdata_chunk
        fft_data_time[n,:] = testdata_time[int(idx[n][0]):int(idx[n][1])]
        fft_data_hr[n,:] = testdata_hr[int(idx[n][0]):int(idx[n][1])]

        # For each data chunk, calaculate the lag_correlations
        nlag = 13 # # of lag-steps to do
        lag_corr = np.zeros(nlag)
        # Use longest data range to computer lag correlations
        for lag in np.arange(0,nlag):
            lag_corr[lag] = np.corrcoef(testdata[int(idx[n][0]):int(idx[n][1]-lag)],testdata[int(lag+idx[n][0]):int(idx[n][1])])[0,1]
        
        #Find time step (Te) where lag correlation reduces to 1/e (0.3678)
        Te[n] = scipy.interpolate.interp1d(lag_corr[:],np.arange(0,nlag),fill_value='extrapolate')(1/np.exp(1)); 
        lagTe[n] = np.exp(-1./Te[n]); # this value is equal to a in ATS655 notes
        lag1[n] = lag_corr[1]

#            fig = plt.figure(figsize=[8,4]);
#            plt.plot(fft_data[n,:],label='Data')
#            plt.plot(fft_data[n,:]*Hanning,label='Data with Hanning')
#            plt.legend()
        
        # Use discrete fast fourier transform to find powers in frequency space 
        temp_arr = np.fft.fft(fft_data[n,:]*Hanning)#, norm = "ortho") #Norm = "ortho" normalizes by 1/np.sqrt(N)
        powers[n,:] = np.abs(temp_arr) # take modulus of data, since fft produces complex values
        power_mod[n,:] = np.abs(powers[n,:]/sizeChunks) # normalize by 1/N    
        power_var[n,:] = power_mod[n,1:int(sizeChunks/2+1)] # Only take indicies 1:N/2, index[0] is zero frequency term (sum of signal / straight horizontal line (not zero size chunk bias many not be equal to the entire dataset)) 
        power_var[n,:] = np.power(power_var[n,:],2.0) # To get power spectrum, need to raise these values (Amplitude Spectrum) to a power of 2
        power_var[n,:] = 2*power_var[n,:] # Multiply by 2 to account for both sides of the spectrum (i.e. the negative frequencies, that were not accounted for in the FFT)
        power_var[n,:] = power_var[n,:]/np.sum(power_var[n,:])
 
        chunk_datetime.append(testdata_time[int(idx[n][0]):int(idx[n][1])])

        test_sum = np.zeros(sizeChunks)
        # Loop through all frequencies and eliminate the power for all other frequencies and recreate the data for each individual frequency
        for f_cnt in np.arange(0,int(sizeChunks)): 
            
            temp_arr_rem4 = copy.deepcopy(temp_arr) # Copy fft data
            
            # Remove all but one frequency at a time            
            temp_arr_rem4[0:f_cnt] = 0; temp_arr_rem4[f_cnt+1:] = 0
            # Convert data back to real data using inverse FFT for each frequency            
            i_data = np.fft.ifft(temp_arr_rem4)
#            print(np.real(i_data))
            test_sum = test_sum + np.real(i_data)    

            # For all frequencies that are less than 1/2 size + 1
            # Find where the peak value is in one period within the dataset
            if f_cnt < int(sizeChunks/2)+1:
#                print(str(per_days[f_cnt-1])+' day cycle')
                one_per = per_hrs[f_cnt-1]/dt # 
            
                ifft_data_re[f_cnt,n,:] = np.real(i_data)
                ifft_data_im[f_cnt,n,:] = np.imag(i_data)
#                plt.plot(np.real(i_data))

                if np.sum(np.isnan(i_data))>0:
                    pk_time[0,f_cnt,n] = np.nan
                    pk_time[1,f_cnt,n] = np.nan
                else:
                    id_imax = np.real(i_data[0:int(one_per)]).argmax()
                    id_imin = np.real(i_data[0:int(one_per)]).argmin()
                    pk_time[0,f_cnt,n] = chunk_datetime[n][id_imax].hour
                    pk_time[1,f_cnt,n] = chunk_datetime[n][id_imin].hour


       
#        # Loop through all frequencies and plot only certain freqeuncies
#        for f_cnt in np.arange(0,int(sizeChunks)): #Start at one since the first frequency is the zero frequency term
#            if f_cnt == 7:
#                print(str(per_days[f_cnt-1])+' day cycle')
#                if np.sum(~np.isnan(ifft_data_re[f_cnt,n,:])) > 0:
#                    import matplotlib.dates as mdates
#
#            
#                    # Add a second frequency
#                    f_cnt2 = 14
#                    print(str(per_days[f_cnt2-1])+' day cycle')
#                    temp_arr_rem4_2 = copy.deepcopy(temp_arr)
#                    temp_arr_rem4_2[0:f_cnt2] = 0
#                    temp_arr_rem4_2[f_cnt2+1:] = 0
#                    i_data_2 = np.fft.ifft(temp_arr_rem4_2)
#                    ifft_data_re_2 = 2*np.real(i_data_2)
#                    
#                    # Add a third frequency
#                    f_cnt3 = 4
#                    print(str(per_days[f_cnt3-1])+' day cycle')
#                    temp_arr_rem4_3 = copy.deepcopy(temp_arr)
#                    temp_arr_rem4_3[0:f_cnt3] = 0
#                    temp_arr_rem4_3[f_cnt3+1:] = 0
#                    i_data_3 = np.fft.ifft(temp_arr_rem4_3)
#                    ifft_data_re_3 = 2*np.real(i_data_3)
#                
#                
#                    fig,ax = plt.subplots(2,1,figsize=(15,8))
#                    
#                    ax[0].plot(chunk_datetime[n],2*ifft_data_re[f_cnt,n,:],'-ob',label=str(per_days[f_cnt-1])+' day cycle')
#                    ax[0].plot(chunk_datetime[n],ifft_data_re_2,'-xg',label=str(per_days[f_cnt2-1])+' day cycle')
#                    ax[0].plot(chunk_datetime[n],ifft_data_re_3,'-dy',label=str(per_days[f_cnt3-1])+' day cycle')
#                    ax[0].plot(chunk_datetime[n],ifft_data_re[f_cnt,n,:]+ifft_data_re_2,'-k',lw=3)
#                    ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%m%d|%H"))
#                    ax[0].set_xlabel('Time')
#                    ax[0].set_ylabel('Real Part of IFFT')
#                    ax[0].set_title('f_cnt:'+str(f_cnt)+' | Chunk:'+str(n))
#                    ax[0].grid()
#                    ax[0].legend(fontsize=10)
#                    
#                    ax2 = plt.twinx(ax[0])
#                    ax2.plot(chunk_datetime[n],test_sum,'-k',lw=5,label='Sum of IFFTs')
#                    ax2.plot(chunk_datetime[n],test_sum/Hanning,':k',lw=5,label='Sum of IFFTs')
#                    ax2.plot(chunk_datetime[n],fft_data[n,:]*Hanning,'-r',lw=3,label='Full Data + Hanning')
#                    ax2.plot(chunk_datetime[n],fft_data[n,:],':r',lw=3,label='Full Data')
#                    ax2.legend(fontsize=10)
#                    ax2.set_ylabel('Full Data')
#
#                    
#                    ax[1].plot(per_days,power_var[n,:])
#                    ax[1].set_xlabel('Period')
#                    ax[1].set_ylabel('Normalizer Power')
#                    ax[1].set_title('f_cnt:'+str(f_cnt)+' | Chunk:'+str(n))
#                    ax[1].set_xscale('log')
#                    ax[1].set_ylim([0,0.3])
#                
##                        print(power_var[n,:])
#                    plt.tight_layout()
##                        #                    plt.savefig('IFFT_JJA_2011_totalmass_f'+str(f_cnt)+'_c'+str(n)+'.png')
##                        #                    plt.close(fig)
##
##        fig = plt.figure(figsize=(10,6))
##        plt.plot(chunk_datetime[n],fft_data[n,:]*Hanning,'-k',label='Original Data')
##        plt.plot(chunk_datetime[n],test_sum,'--y',label='Sum of Frequencies')
##        plt.legend()

    cnt_Nchunks_yy = np.count_nonzero(~np.isnan(power_var[:,1])) # count number if non-NAN chunks per year
    return power_var,freq,per_days,pk_time[0,:,:],pk_time[1,:,:],ifft_data_re,ifft_data_im, lagTe, lag1, lag1_pd, chunk_datetime

def red_noise(lag_value,freq_in):        
    red_lag = (1-np.power(lag_value,2)) / (1 - (2*lag_value*np.cos(freq_in*2*np.pi)) + np.power(lag_value,2));
    red_lag_n = red_lag/np.sum(red_lag);
    return red_lag, red_lag_n



#-------------------------------------------------------------------------------------------------------------------------


##############################################################################
##############################################################################
#####
#####           Read-In Data
#####
##############################################################################
##############################################################################

savepath = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/'
#savepath = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/H/'
df_120M = pickle.load(open(savepath+"df_120M_FnvFinal_D2E4_var95_LT2wk18_CPCQC.p","rb"))
df_120M = df_120M[0]

# Fill data points that are missing three points in a row
for i in np.arange(0,np.shape(df_120M)[1]):
    cur_values = copy.deepcopy(df_120M.iloc[:,i].values)
    idxs = []
    for j in np.arange(1,len(cur_values)-4):
        if np.isnan(cur_values[j]) and np.isnan(cur_values[j+1]) and np.isnan(cur_values[j+2]) and ~np.isnan(cur_values[j-1]) and ~np.isnan(cur_values[j+3]):
            idxs.append(j)
    new_cur_values = copy.deepcopy(cur_values)
    for j in np.arange(0,len(idxs)):
        slope = (new_cur_values[idxs[j]+3] - new_cur_values[idxs[j]-1])/4
        new_cur_values[idxs[j]] = (new_cur_values[idxs[j]-1]+(1*slope))
        new_cur_values[idxs[j]+1] = (new_cur_values[idxs[j]-1]+(2*slope))
        new_cur_values[idxs[j]+2] = (new_cur_values[idxs[j]-1]+(3*slope))
        if np.abs(new_cur_values[idxs[j]+3]-(new_cur_values[idxs[j]-1]+(4*slope))) > 0.001:
            dgsgsdfgfsdgfdgfdgfsgersdfgsdg
    df_120M.iloc[:,i] = copy.deepcopy(new_cur_values)

# Fill data points that are missing two points in a row
for i in np.arange(0,np.shape(df_120M)[1]):
    cur_values = copy.deepcopy(df_120M.iloc[:,i].values)
    idxs = []
    for j in np.arange(1,len(cur_values)-3):
        if np.isnan(cur_values[j]) and np.isnan(cur_values[j+1]) and ~np.isnan(cur_values[j-1]) and ~np.isnan(cur_values[j+2]):
            idxs.append(j)
    new_cur_values = copy.deepcopy(cur_values)
    for j in np.arange(0,len(idxs)):
        slope = (new_cur_values[idxs[j]+2] - new_cur_values[idxs[j]-1])/3
        new_cur_values[idxs[j]] = (new_cur_values[idxs[j]-1]+(1*slope))
        new_cur_values[idxs[j]+1] = (new_cur_values[idxs[j]-1]+(2*slope))
        if np.abs(new_cur_values[idxs[j]+2]-(new_cur_values[idxs[j]-1]+(3*slope))) > 0.001:
            dgsgsdfgfsdgfdgfdgfsgersdfgsdg
    df_120M.iloc[:,i] = copy.deepcopy(new_cur_values)

# Fill data points that are missing one point in a row
for i in np.arange(0,np.shape(df_120M)[1]):
    cur_values = copy.deepcopy(df_120M.iloc[:,i].values)
    idxs = []
    for j in np.arange(1,len(cur_values)-2):
        if np.isnan(cur_values[j]) and ~np.isnan(cur_values[j-1]) and ~np.isnan(cur_values[j+1]):
            idxs.append(j)
    new_cur_values = copy.deepcopy(cur_values)
    for j in np.arange(0,len(idxs)):
        slope = (new_cur_values[idxs[j]+1] - new_cur_values[idxs[j]-1])/2
        new_cur_values[idxs[j]] = (new_cur_values[idxs[j]-1]+(1*slope))
        if np.abs(new_cur_values[idxs[j]+1]-(new_cur_values[idxs[j]-1]+(2*slope))) > 0.001:
            dgsgsdfgfsdgfdgfdgfsgersdfgsdg
    df_120M.iloc[:,i] = copy.deepcopy(new_cur_values)

# Group data by season
dfs = ['120M']; mo_num = 23;
for i in np.arange(0,len(dfs)):
    exec('df_cur = df_'+dfs[i])
    df_cur_MAM = df_cur[(df_cur.iloc[:,mo_num] == 3) | (df_cur.iloc[:,mo_num] == 4) | (df_cur.iloc[:,mo_num] == 5)]
    df_cur_JJA = df_cur[(df_cur.iloc[:,mo_num] == 6) | (df_cur.iloc[:,mo_num] == 7) | (df_cur.iloc[:,mo_num] == 8)]
    df_cur_SON = df_cur[(df_cur.iloc[:,mo_num] == 9) | (df_cur.iloc[:,mo_num] == 10) | (df_cur.iloc[:,mo_num] == 11)]
    df_cur_D = df_cur[(df_cur.iloc[:,mo_num] == 12)]
    df_cur_JF = df_cur[(df_cur.iloc[:,mo_num] == 1) | (df_cur.iloc[:,mo_num] == 2)]
    exec('df_'+dfs[i]+'_MAM = df_cur_MAM')
    exec('df_'+dfs[i]+'_JJA = df_cur_JJA')
    exec('df_'+dfs[i]+'_SON = df_cur_SON')
    exec('df_'+dfs[i]+'_JF = df_cur_JF')
    exec('df_'+dfs[i]+'_D = df_cur_D')
df_120M_DJF = pd.concat([df_120M_D,df_120M_JF])

# Group seasonal data by year
dfs = ['120M']; yr_num = mo_num - 1
yys = np.arange(2009,2015)
for i in np.arange(0,len(dfs)):
    for j in yys: 
        exec('df_'+dfs[i]+'_'+str(j)+' = df_'+dfs[i]+'[(df_'+dfs[i]+'.iloc[:,'+str(yr_num)+'] == '+str(j)+')]')
        exec('df_'+dfs[i]+'_MAM_'+str(j)+' = df_'+dfs[i]+'_MAM[(df_'+dfs[i]+'_MAM.iloc[:,'+str(yr_num)+'] == '+str(j)+')]')
        exec('df_'+dfs[i]+'_JJA_'+str(j)+' = df_'+dfs[i]+'_JJA[(df_'+dfs[i]+'_JJA.iloc[:,'+str(yr_num)+'] == '+str(j)+')]')
        exec('df_'+dfs[i]+'_SON_'+str(j)+' = df_'+dfs[i]+'_SON[(df_'+dfs[i]+'_SON.iloc[:,'+str(yr_num)+'] == '+str(j)+')]')
        exec('df_'+dfs[i]+'_DJF_'+str(j)+' = pd.concat( [ df_'+dfs[i]+'_D[(df_'+dfs[i]+'_D.iloc[:,'+str(yr_num)+'] == '+str(j-1)+')] , df_'+dfs[i]+'_JF[(df_'+dfs[i]+'_JF.iloc[:,'+str(yr_num)+'] == '+str(j)+')] ])')       

savepath = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/'
##############################################################################
##############################################################################
#####
#####            Specify Inputs
#####
##############################################################################
##############################################################################
casenames = ['120M','120M_MAM','120M_JJA','120M_SON','120M_DJF']; 
colors=['black','blue','crimson','goldenrod','0.75']

varnames = ['intN0']
#varnames = ['intN_v2','intN0','intN1','intN2','intN3']
#varnames = ['intS_v2','intS0','intS1','intS2','intS3']
#varnames = ['intV_v2','intV0','intV1','intV2','intV3']
#varnames = ['intN_v2']

cnames_plt = ['MAM','JJA','SON','DJF']
dt = 2
num_days = 7
overlap = 0.5
days_conv = 24 #convert hourly data to daily data

###########
#### CHANGE BELOW TO CREATE DIFFERENT PLOTS
###########
plot_combined = 0; 
if plot_combined == 1:
    varnames = ['intN0','intN1','intN2','intN3']
#    varnames = ['intS0','intS1','intS2','intS3']
#    varnames = ['intV0','intV1','intV2','intV3']
    plt.rcParams.update({'font.size': 14})
    fig, axs = plt.subplots(5,4,figsize=(18,14))
elif plot_combined == 2:
#    varnames = ['intN0','intN1','intN2','intN3']
#    varnames = ['intS0','intS1','intS2','intS3']
#    varnames = ['intV0','intV1','intV2','intV3']
    plt.rcParams.update({'font.size': 14})
    fig, axs = plt.subplots(4,1,figsize=(7,12))
   
##############################################################################
##############################################################################
#####
#####            Combine data from different years and cases
#####     [ Create red noise spectra, significance testing, and plotting]
##############################################################################
##############################################################################

for v in np.arange(0,len(varnames)):
    varname = varnames[v]

    addsave = 'x7day_xxxxPUB_'+varname+'_'+str(num_days)+'day_'+str(dt*60)+'M_'+str(overlap)+'overlap_ylim'

    # Initiliaze Arrays 
    pwr_s_comb = []; maxt_s_comb = []; mint_s_comb = []
    ifft_re_s_comb = []; ifft_im_s_comb = []
    lagTe_s_comb = []; lag1_s_comb = []
    lag1_pd_s_comb = []
    chunktime_s_comb = []; Te_s_comb = []

    for case_num in np.arange(0,len(casenames)):
        cur_case = casenames[case_num]   
        cnt_yy = 0
        for yyyy in np.arange(2009,2014):   
            exec('cur_data = df_'+cur_case+'_'+str(yyyy)+'.'+varname)
            [pwr_now,frq,pr_dys,max_times,min_times,ifft_re,ifft_im,lagTe,lag1,lag1_pd,chunk_times] = power_spectrum(cur_data,dt,num_days,overlap) # Data, dt (hours), chunk size (days), overlap (fraction)
            
            if np.size(pwr_now) > 0:
            
                Te = 1 / np.log(lagTe)
                if cnt_yy == 0:
                    pwr_comb = pwr_now; maxt_comb = max_times; mint_comb = min_times
                    ifft_re_comb = ifft_re; ifft_im_comb = ifft_im; 
                    lagTe_comb = lagTe; lag1_comb = lag1; lag1_pd_comb = lag1_pd; 
                    ctime_comb = chunk_times; Te_comb = Te;
                else:
                    pwr_comb = np.vstack((pwr_comb,pwr_now))
                    maxt_comb = np.hstack((maxt_comb,max_times))
                    mint_comb = np.hstack((mint_comb,min_times))
                    ifft_re_comb = np.hstack((ifft_re_comb,ifft_re));
                    ifft_im_comb = np.hstack((ifft_im_comb,ifft_im));
                    lagTe_comb = np.hstack((lagTe_comb,lagTe)); lag1_comb = np.hstack((lag1_comb,lag1)); lag1_pd_comb = np.hstack((lag1_pd_comb,lag1_pd)); 
                    ctime_comb = np.vstack((ctime_comb,chunk_times)); Te_comb = np.hstack((Te_comb,Te));
                
                cnt_yy = cnt_yy + 1
    
        pwr_s_comb.append(pwr_comb) # power associated with frequencies
        maxt_s_comb.append(maxt_comb) # peak in cycle 
        mint_s_comb.append(mint_comb) # peak in cycle
        ifft_re_s_comb.append(ifft_re_comb )
        ifft_im_s_comb.append(ifft_im_comb )
        lagTe_s_comb.append(lagTe_comb )
        lag1_s_comb.append(lag1_comb )
        lag1_pd_s_comb.append(lag1_pd_comb )
        chunktime_s_comb.append(ctime_comb )
        Te_s_comb.append(Te_comb)

    ##############################################################################
    ##############################################################################
    #####
    #####            Calculate red noise spectrum and significance
    #####
    ##############################################################################
    ##############################################################################
    from scipy.stats import f, chi2
    red_spec = [];
    F99red1 = []; F99red2 = []; F99red3 = []
    lag_value = np.zeros(len(casenames))
    DOF = np.zeros(len(casenames))
    red_noise_norm = np.zeros((5,len(frq)))
    for case_num in np.arange(0,len(casenames)):
    
        num_chunks = np.sum(~np.isnan(lag1_s_comb[case_num]))
        #print(num_chunks)       
        cur_case = casenames[case_num]
    
        # Uses either e-folding or lag-1 autocorrelations
        lag_value[case_num] = np.nanmedian(lag1_s_comb[case_num])  # Lag-1 autocor
        #lag_value[case_num] = np.nanmedian(lagTe_s_comb[case_num]) # e-folding time
       
        # Calculate red noise spectrum, using median lag value
        red_spec.append(red_noise(lag_value[case_num],frq))
    
        # Calculate red noise spectrum from each individual lag value, and calculate mean or median of all the red noise spectra
        all_red_noise_norm = []
        for k in np.arange(0,len(lag1_s_comb[case_num])):
            if ~np.isnan(lag1_s_comb[case_num][k]):
                temp = red_noise(lag1_s_comb[case_num][k],frq)
                all_red_noise_norm.append(temp[1])
        red_noise_norm[case_num,:] = np.nanmedian(np.asarray(all_red_noise_norm),axis=0)
        
        sig_level = 0.99
        #sig_level = np.power(0.99,1.2/(sizeChunks/2)) # 1.2 to count for smoothing, still want sig_level but taking into account a posteori statistics
        DOF[case_num] = 2.8*num_chunks #Welch 1967 (effective Dof)
        F99red1.append(f.ppf(sig_level,DOF[case_num],1000) * red_spec[case_num][1]) #Uses mean/median value, 1 in red_spec pulls out normalized values
        F99red2.append(f.ppf(sig_level,DOF[case_num],1000) * red_noise_norm[case_num,:]) #Uses mean/median red noise spectrum, calculatd from individual red noise 
        #F99red3.append(chi2.ppf(sig_level,DOF[case_num]) * red_noise_norm[case_num,:] / DOF[case_num]) # red noise spectrum averaged for individual red noise 
    
    
    ####### Test red noise spectrum function
    #for i in np.arange(0,5):
    #    fig = plt.figure(figsize=(8,5))
    #    plt.plot(frq,red_spec[i][1],label=str(i)+' Mean Lag1')
    #    plt.plot(frq,red_noise_norm[i,:],label=str(i)+' Mean Spectra')
    #    plt.legend()
    #plt.xlabel('Frequency')
    #plt.ylabel('Normalized Power')
    #
    
    ####### Test red noise spectrum function
    #fig = plt.figure(figsize=(8,5))
    #for i in np.arange(0.5,0.91,0.1):
    #    rnnow = red_noise(i,frq)
    #    plt.plot(frq,rnnow[1],label=str(i))
    #plt.legend()
    #plt.xlabel('Frequency')
    #plt.ylabel('Normalized Power')
    
    ##########################
    ####### PLOT
    ##########################

    if plot_combined == 0:       

        from cycler import cycler
        plt.rcParams['axes.prop_cycle'] = cycler(color=['black','blue','crimson','goldenrod','0.75'])

        ########################################################################
        ########################################################################
        ######### Plot power spectrum for all the cases on ONE FIGURE
        ########################################################################
        ########################################################################

        plt.rcParams.update({'font.size': 15})
        fig, axs = plt.subplots(1,1,figsize=(7,5))
        axs.plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[0],axis=0),'-o',lw=4,label='All ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[0][:,0])))+')')
        axs.plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[1],axis=0),'-',lw=2,label='MAM ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[1][:,0])))+')')
        axs.plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[2],axis=0),'-',lw=2,label='JJA ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[2][:,0])))+')')
        axs.plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[3],axis=0),'-',lw=2,label='SON ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[3][:,0])))+')')
        axs.plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[4],axis=0),'-',lw=2,label='DJF ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[4][:,0])))+')')
        
        #axs.plot(pr_dys*days_conv,red_spec[0][1],':',lw=2,label='All Red Noise')
        #axs.plot(pr_dys*days_conv,red_spec[1][1],':',lw=2)
        #axs.plot(pr_dys*days_conv,red_spec[2][1],':',lw=2)
        #axs.plot(pr_dys*days_conv,red_spec[3][1],':',lw=2)
        #axs.plot(pr_dys*days_conv,red_spec[4][1],':',lw=2)
        
        axs.plot(pr_dys*days_conv,red_noise_norm[0,:],':',lw=2,label='All Red Noise')
        axs.plot(pr_dys*days_conv,red_noise_norm[1,:],':',lw=2)
        axs.plot(pr_dys*days_conv,red_noise_norm[2,:],':',lw=2)
        axs.plot(pr_dys*days_conv,red_noise_norm[3,:],':',lw=2)
        axs.plot(pr_dys*days_conv,red_noise_norm[4,:],':',lw=2)
        
        axs.plot(pr_dys*days_conv,F99red2[0],':',lw=1,label='99% Confidence Level')
        axs.plot(pr_dys*days_conv,F99red2[1],':',lw=1)
        axs.plot(pr_dys*days_conv,F99red2[2],':',lw=1)
        axs.plot(pr_dys*days_conv,F99red2[3],':',lw=1)
        axs.plot(pr_dys*days_conv,F99red2[4],':',lw=1)
        
        #axs.set_yscale('log')
        axs.set_xscale('log')
        #axs.set_xticks([0.25,0.5,1,3]*24,minor = False)
        #axs.set_xticks([6,12,24,72,168],minor = False)
        axs.grid(which = 'major')
        
        axs.xaxis.set_major_formatter(FormatStrFormatter('%i'))
        axs.set_xlabel('Period (Hours)')
        axs.set_ylabel('Normalized Power')
        axs.set_ylim([0,0.2])
        #axs.set_xlim([dt*2/24*days_conv,num_days*days_conv/2])
        
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(savepath+addsave+'Power_Spectrum_'+varname+'_Narrow.png')
        plt.savefig(savepath+addsave+'Power_Spectrum_'+varname+'_Narrow.eps')
        plt.close(fig)
    
        ########################################################################
        ########################################################################
        ######### Plot power spectrum for all the cases on SEPARATE Figures
        ########################################################################
        ########################################################################
        
        plt.rcParams.update({'font.size': 15})
        fig, axs = plt.subplots(5,1,figsize=(6,13))
        colors=['black','blue','crimson','goldenrod','0.75']
        texts = ['a) All','b) MAM','c) JJA','d) SON','e) DJF']
        for k in np.arange(0,5):
            #axs[k].plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[k],axis=0),'o',color=colors[k],label='All ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[k][:,0])))+')')
            #axs[k].plot(pr_dys*days_conv,red_spec[k][1],':',color=colors[k],lw=4)
            #axs[k].plot(pr_dys*days_conv,F99red[k],'-',color=colors[k],lw=3)
            axs[k].plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[k],axis=0),'o',color=colors[k],label='All ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[k][:,0])))+')')

        #    axs[k].plot(pr_dys*days_conv,red_spec[k][1],':',color=colors[k],lw=4)
        #    axs[k].plot(pr_dys*days_conv,F99red1[k],'-',color=colors[k],lw=3)
            axs[k].plot(pr_dys*days_conv,red_noise_norm[k,:],':',color=colors[k],lw=4)
            axs[k].plot(pr_dys*days_conv,F99red2[k],'-',color=colors[k],lw=3)
            axs[k].set_xlim([2,14])
            axs[k].set_xscale('log')
            axs[k].grid(which = 'major')
            
            axs[k].xaxis.set_major_formatter(FormatStrFormatter('%i'))

            #SMPS       
            axs[k].set_xticks([4,6,12,24,48,84])
            axs[k].set_xticklabels([4,6,12,24,48,84])   
            axs[k].set_ylim([0.001,0.2])
            axs[k].set_yticks([0,0.05,0.1,0.15,0.2])   
            axs[k].set_yticklabels(['0','0.05','0.10','0.15','0.20'])   

#            #CPC
#            axs[k].set_xticks([2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
#            axs[k].set_xticklabels([2,3,4,5,6,7,8,'',10,'',12,'',14])   
#            axs[k].set_ylim([0.001,0.22])
#            axs[k].set_yticks([0,0.05,0.1,0.15,0.2])   
#            axs[k].set_yticklabels(['0','0.05','0.10','0.15','0.20'])   
        
            axs[k].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            axs[k].set_ylabel('Normalized Power')   
            axs[k].set_xlim([dt*2/24*days_conv * 0.97,num_days*days_conv/2 * 1.03])
            #axs[k].set_xlim([5.5,num_days*days_conv/2 * 1.03])
            
        #    axs[k].text(2.3,0.12,texts[k],fontsize = 18)
        
            #SMPS
            axs[k].text(4.4,0.165,texts[k]+' ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[k][:,0])))+')',fontsize = 15)
#            #CPC
#            axs[k].text(2.15,0.16,texts[k]+' ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[k][:,0])))+')',fontsize = 18)
        
            if k == 4:
                axs[k].set_xlabel('Period (Hours)')
#                axs[k].set_xlabel('Period (Days)')
        
        #plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(savepath+addsave+'Power_Spectrum_'+varname+'_Separate.png')
        plt.savefig(savepath+addsave+'Power_Spectrum_'+varname+'_Separate.eps')
        plt.close(fig)
    
        ########################################################################
        ########################################################################
        ######### Plot histogram of peak times for each cycle
        ########################################################################
        ########################################################################
        if cur_case.__contains__('1440'):

            find_freq_4d = np.abs(pr_dys-4.0).argmin() + 1
            #print(str(pr_dys[find_freq_24-1])+' day cycle')
            hours_for_analysis = pr_dys[find_freq_4d-1];
            
            #pwr_thresh = 0.05
            cnt_plt = 0

            for i in np.arange(1,5):            
                
                # Choose power thresholdhold to be greater than the 99th percentile above red noise for each individual Season
                pwr_thresh = copy.deepcopy(F99red2[i][find_freq_4d-1])
                pwr_thresh = copy.deepcopy(red_noise_norm[i,find_freq_4d-1])
                #print(pwr_thresh)
                
                temp_times = chunktime_s_comb[i]
                temp_data = copy.deepcopy(maxt_s_comb[i][find_freq_4d,:])
                temp_pwr_data = copy.deepcopy(pwr_s_comb[i][:,find_freq_4d])
                for ii in np.arange(0,len(temp_data)):
                    if temp_pwr_data[ii] < pwr_thresh or np.isnan(temp_pwr_data[ii]):
                        temp_data[ii] = np.nan
            
                # Save times where there is a lot of power associated with the 4d cycle
                tp = []; tp_data = []; tp_4ddata = []; tp_pwr = []
                for k in np.arange(0,len(temp_data)):
                    if ~np.isnan(temp_data[k]):
                        tp.append(temp_times[k,:]) # times powerful
                        tp_4ddata.append(ifft_re_s_comb[i][find_freq_4d,k,:])
                        tp_data.append(np.nansum(ifft_re_s_comb[i][:,k,:],axis=0))
                        tp_pwr.append(temp_pwr_data[k])
    
    
                # Save times where there is no power associated with the diurnal cycle              
                tnp = []; tnp_data = []; tnp_4ddata = []; tnp_pwr = []
                for k in np.arange(0,len(temp_data)):
                    if temp_pwr_data[k] < red_noise_norm[i,find_freq_4d]:
                        tnp.append(temp_times[k,:]) # times powerful
                        tnp_4ddata.append(ifft_re_s_comb[i][find_freq_4d,k,:])
                        tnp_data.append(np.nansum(ifft_re_s_comb[i][:,k,:],axis=0))
                        tnp_pwr.append(temp_pwr_data[k])
    
                cnt_vals = np.sum(np.logical_not(np.isnan(temp_data)))
                temp_data = temp_data[np.logical_not(np.isnan(temp_data))]
    
                #
                import pickle
                pickle.dump((tp,tp_data,tp_4ddata,tp_pwr,tnp,tnp_data,tnp_4ddata,tnp_pwr),open(savepath+"times_power_"+cnames_plt[cnt_plt]+"_"+varname+"_4day.p","wb"))
                cnt_plt = cnt_plt + 1
                
        
        if cur_case.__contains__('120'):
            #############################################################
            #### 24 hour cycle
            #############################################################
            find_freq_24 = np.abs(pr_dys-1.0).argmin() + 1
            #print(str(pr_dys[find_freq_24-1])+' day cycle')
            hours_for_analysis = pr_dys[find_freq_24-1]*24;
            
            #pwr_thresh = 0.05
            
            x_utc = np.arange(0,25,2) # local daylight time
            x_utc = np.insert(x_utc,13,24.1)
            x_cdt = np.arange(0,25,2)-5 # utc time
            for i in np.arange(0,len(x_cdt)):
                if x_cdt[i] < 0:
                    x_cdt[i] = x_cdt[i]+24            
            
            plt.rcParams.update({'font.size': 15})
            fig, axs = plt.subplots(1,1,figsize=(9,5))
            cnt_plt = 0
            max_val = 0
            next(axs._get_lines.prop_cycler) # Skip the black line for ALL for this plot.
            for i in np.arange(1,5):
                
                # Choose power thresholdhold to be greater than the 99th percentile above red noise for each individual Season
                pwr_thresh = copy.deepcopy(F99red2[i][find_freq_24-1])
                pwr_thresh = copy.deepcopy(red_noise_norm[i,find_freq_24-1])
                #print(pwr_thresh)
                
                temp_times = chunktime_s_comb[i]
                temp_data = copy.deepcopy(maxt_s_comb[i][find_freq_24,:])
                temp_pwr_data = copy.deepcopy(pwr_s_comb[i][:,find_freq_24])
                for ii in np.arange(0,len(temp_data)):
                    if temp_pwr_data[ii] < pwr_thresh or np.isnan(temp_pwr_data[ii]):
                        temp_data[ii] = np.nan
            
                # Save times where there is a lot of power associated with the diurnal cycle
                tp = []; tp_data = []; tp_24data = []; tp_pwr = []
                for k in np.arange(0,len(temp_data)):
                    if ~np.isnan(temp_data[k]):
                        tp.append(temp_times[k,:]) # times powerful
                        tp_24data.append(ifft_re_s_comb[i][find_freq_24,k,:])
                        tp_data.append(np.nansum(ifft_re_s_comb[i][:,k,:],axis=0))
                        tp_pwr.append(temp_pwr_data[k])
    
    
                # Save times where there is no power associated with the diurnal cycle              
                tnp = []; tnp_data = []; tnp_24data = []; tnp_pwr = []
                for k in np.arange(0,len(temp_data)):
                    if temp_pwr_data[k] < red_noise_norm[i,find_freq_24]:
                        tnp.append(temp_times[k,:]) # times powerful
                        tnp_24data.append(ifft_re_s_comb[i][find_freq_24,k,:])
                        tnp_data.append(np.nansum(ifft_re_s_comb[i][:,k,:],axis=0))
                        tnp_pwr.append(temp_pwr_data[k])
    
                cnt_vals = np.sum(np.logical_not(np.isnan(temp_data)))
                temp_data = temp_data[np.logical_not(np.isnan(temp_data))]
    
                #
                import pickle
                pickle.dump((tp,tp_data,tp_24data,tp_pwr,tnp,tnp_data,tnp_24data,tnp_pwr),open(savepath+"times_power_"+cnames_plt[cnt_plt]+"_"+varname+"_7day.p","wb"))
            
                bins = np.arange(0,25,2) 
                a1,a2 = np.histogram(temp_data,bins=bins)
                a1 = np.insert(a1,0,0)
                a1 = np.insert(a1,13,0)
                a1_n = a1/np.sum(a1)
            
                axs.plot(x_utc,a1_n,ls='steps',lw=5-(0.5*i),label=cnames_plt[cnt_plt]+'('+str(cnt_vals)+')')
                axs.set_xticks(x_utc)
            
                cnt_plt = cnt_plt + 1
            
                max_val = np.max([np.max(a1_n),max_val])
            
            axs.legend(fontsize=13,loc='best')
            axs.set_ylabel('Normalized Frequency')
            axs.set_xlabel('Hour of Day (UTC Time)')

            plt.rcParams['xtick.labelsize']=15
            ax2 = axs.twiny()
            ax2.set_xticks(x_utc)
            ax2.set_xticklabels(x_cdt,fontsize=15)
            ax2.xaxis.set_ticks_position('bottom')
            ax2.xaxis.set_label_position('bottom')
            ax2.spines['bottom'].set_position(('outward',50))
            ax2.set_xlabel('Hour of Day (Central Daylight Time)',fontsize=15)
            ax2.set_xlim(axs.get_xlim())

            plt.tight_layout()
            plt.ylim([0,max_val+0.05])
            plt.savefig(savepath+addsave+'Peak_'+str(pr_dys[find_freq_24-1])+'day_Cycle_Histogram_'+varname+'_'+str(round(pwr_thresh,3))+'pwrthresh.png')
            plt.savefig(savepath+addsave+'Peak_'+str(pr_dys[find_freq_24-1])+'day_Cycle_Histogram_'+varname+'_'+str(round(pwr_thresh,3))+'pwrthresh.eps')
            plt.close(fig)
            
            #############################################################
            #### 12 hour cycle
            #############################################################
            find_freq = np.abs(pr_dys-0.5).argmin() + 1
            #print(str(pr_dys[find_freq-1])+' day cycle')
            hours_for_analysis = pr_dys[find_freq-1]*24;
    #        
            plt.rcParams.update({'font.size': 15})
            fig, axs = plt.subplots(1,1,figsize=(9,5))
            cnt_plt = 0
            max_val = 0
            next(axs._get_lines.prop_cycler) # Skip the black line for all for this plot.
            for i in np.arange(1,5):#5):
                
                # Choose power thresholdhold to be greater than the 99th percentile above red noise for each individual Season
                pwr_thresh = copy.deepcopy(F99red2[i][find_freq-1])
                pwr_thresh = copy.deepcopy(red_noise_norm[i,find_freq-1])
                #print(pwr_thresh)
                
                temp_times = chunktime_s_comb[i]
                temp_data = copy.deepcopy(maxt_s_comb[i][find_freq,:])
                temp_pwr_data = copy.deepcopy(pwr_s_comb[i][:,find_freq])
                for ii in np.arange(0,len(temp_data)):
                    if temp_pwr_data[ii] < pwr_thresh or np.isnan(temp_pwr_data[ii]):
                        temp_data[ii] = np.nan
                                    
    #            # Save times where there is a lot of power associated with the diurnal cycle
    #            times_powerful_12 = []
    #            ifft_re_1212 = []
    #            ifft_re_1224 = []
    #            ifft_im_1212 = []
    #            ifft_im_1224 = []
    #            for k in np.arange(0,len(temp_data)):
    #                if ~np.isnan(temp_data[k]):
    #                    times_powerful_12.append(temp_times[k,:])
    #
    #                    # Plot inverse ffts for each time period
    #                    cur_ifft_im = ifft_im_s_comb[i][find_freq,k,:]
    #                    cur_ifft_re = ifft_re_s_comb[i][find_freq,k,:]
    #
    #                    ifft_im_1212.append(cur_ifft_im)
    #                    ifft_re_1212.append(cur_ifft_re)
    #
    #                    # Plot inverse ffts for each time period
    #                    cur_ifft_im_24 = ifft_im_s_comb[i][find_freq_24,k,:]
    #                    cur_ifft_re_24 = ifft_re_s_comb[i][find_freq_24,k,:]
    #
    ##                    ifft_im_1224 = np.append(ifft_im_1224,cur_ifft_im_24,axis=i-1)
    ##                    ifft_re_1224 = np.append(ifft_re_1224,cur_ifft_re_24,axis=i-1)
    #
    #                    ifft_im_1224 = np.hstack((ifft_im_1224,cur_ifft_im_24))
    #                    ifft_re_1224 = np.hstack((ifft_re_1224,cur_ifft_re_24))
    #
    #
    ##                    import matplotlib.dates as mdates
    #
    #                    fig,ax = plt.subplots(1,1,figsize=(20,5))
    #                    #plt.plot(temp_times[k,:],2*cur_ifft_im + 2*cur_ifft_im_24,label='Im')
    #                    plt.plot(temp_times[k,:],2*cur_ifft_re,label='12 Hour Cycle')
    #                    plt.plot(temp_times[k,:],2*cur_ifft_re_24,label='24 Hour Cycle')
    #                    
    #                    plt.plot(temp_times[k,:],2*cur_ifft_re+2*cur_ifft_re_24,label='Combined')
    #                    plt.legend()
    #                    plt.title(cnames_plt[cnt_plt]+' '+varname+'_'+str(temp_data[k])+' '+str(temp_pwr_data[k]))
    #                    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
    
    #        dfggsfg
                # Save times where there is no power associated with the diurnal cycle              
    #            times_notpowerful = []
    #            for k in np.arange(0,len(temp_data)):
    #                if temp_pwr_data[k] < red_noise_norm[i,find_freq]:
    #                    times_notpowerful.append(temp_times[k,:])
    #
    #            cnt_vals = np.sum(np.logical_not(np.isnan(temp_data)))
    #            temp_data = temp_data[np.logical_not(np.isnan(temp_data))]
    #
    #            #
    #            pathname = '/Volumes/Data/Act_Spectrum/'
    #            import pickle
    #            pickle.dump((times_powerful,temp_data,times_notpowerful),open(pathname+"times_power_12hr_"+cnames_plt[cnt_plt]+"_"+varname+".p","wb"))
    
                cnt_vals = np.sum(np.logical_not(np.isnan(temp_data)))
                temp_data = temp_data[np.logical_not(np.isnan(temp_data))]
            
                bins = np.arange(0,25,2) 
                a1,a2 = np.histogram(temp_data,bins=bins)
                a1 = np.insert(a1,0,0)
                a1 = np.insert(a1,13,0)
                a1_n = a1/np.sum(a1)
            
                a1_temp = (a1_n[1:7]+a1_n[7:13])
                a1_n[1:7] = a1_temp
                a1_n[7:13] = a1_temp
            
                axs.plot(x_utc,a1_n,ls='steps',lw=5-(0.5*i),label=cnames_plt[cnt_plt]+'('+str(cnt_vals)+')')
                axs.set_xticks(x_utc)
            
                cnt_plt = cnt_plt + 1
            
                max_val = np.max([np.max(a1_n),max_val])
                        
            axs.legend(fontsize=13)
            axs.set_ylabel('Normalized Frequency')
            axs.set_xlabel('Hour of Day (UTC Time)')

            plt.rcParams['xtick.labelsize']=15
            ax2 = axs.twiny()
            ax2.set_xticks(x_utc)
            ax2.set_xticklabels(x_cdt,fontsize=15)
            ax2.xaxis.set_ticks_position('bottom')
            ax2.xaxis.set_label_position('bottom')
            ax2.spines['bottom'].set_position(('outward',50))
            ax2.set_xlabel('Hour of Day (Central Daylight Time)',fontsize=15)
            ax2.set_xlim(axs.get_xlim())

            plt.tight_layout()
            plt.ylim([0,max_val+0.05])
            plt.savefig(savepath+addsave+'Peak_'+str(pr_dys[find_freq-1])+'day_Cycle_Histogram_'+varname+'_'+str(round(pwr_thresh,3))+'pwrthresh.png')
            plt.savefig(savepath+addsave+'Peak_'+str(pr_dys[find_freq-1])+'day_Cycle_Histogram_'+varname+'_'+str(round(pwr_thresh,3))+'pwrthresh.eps')
            plt.close(fig)
        
        
    elif plot_combined == 1:       
        
        from cycler import cycler
        plt.rcParams['axes.prop_cycle'] = cycler(color=['black','blue','crimson','goldenrod','0.75'])

        plt.rcParams.update({'font.size': 15})

        if v == 0 :
            texts = ['a) All','b) MAM','c) JJA','d) SON','e) DJF']
        elif v == 1 :
            texts = ['f) All','g) MAM','h) JJA','i) SON','j) DJF']
        elif v == 2 :
            texts = ['k) All','l) MAM','m) JJA','n) SON','o) DJF']
        elif v == 3 :
            texts = ['p) All','q) MAM','r) JJA','s) SON','t) DJF']

        for k in np.arange(0,5):
            axs[k,v].plot(pr_dys*days_conv,np.nanmean(pwr_s_comb[k],axis=0),'o',color=colors[k],label='All ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[k][:,0])))+')')
            axs[k,v].plot(pr_dys*days_conv,red_noise_norm[k,:],':',color=colors[k],lw=4)
            axs[k,v].plot(pr_dys*days_conv,F99red2[k],'-',color=colors[k],lw=3)
            axs[k,v].set_xlim([2,14])
            axs[k,v].set_xscale('log')
            axs[k,v].grid(which = 'major')       
            axs[k,v].xaxis.set_major_formatter(FormatStrFormatter('%i'))
            
            axs[k,v].set_xticks([4,6,12,24,48,84])
            axs[k,v].set_xticklabels([4,6,12,24,48,84],fontsize=15)   
            axs[k,v].set_ylim([0.001,0.2])
            axs[k,v].set_yticks([0,0.05,0.1,0.15,0.2])
            axs[k,v].set_yticklabels(['0','0.05','0.10','0.15','0.20'],fontsize=15)   
            
            axs[k,v].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            axs[k,v].set_xlim([dt*2/24 * days_conv * 0.97,num_days*days_conv/2 * 1.03])
            #axs[k,v].set_xlim([5.5,num_days*days_conv/2 * 1.03])

            axs[k,v].text(4.4,0.165,texts[k]+' ('+str(np.count_nonzero(~np.isnan(pwr_s_comb[k][:,0])))+')',fontsize = 15)

            if k == 0 and v == 0:
                axs[k,v].set_title('7-30nm',fontsize=16)

            if k == 0 and v == 1:
                axs[k,v].set_title('30-140nm',fontsize=16)

            if k == 0 and v == 2:
                axs[k,v].set_title('140-800nm',fontsize=16)

            if k == 0 and v == 3:
                axs[k,v].set_title('800nm+',fontsize=16)
                
            if v == 0:    
                axs[k,v].set_ylabel('Normalized Power',fontsize=16)   
            if k == 4:
                axs[k,v].set_xlabel('Period (Hours)',fontsize=16)
                
                
    elif plot_combined == 2:       
        
        from cycler import cycler
        plt.rcParams['axes.prop_cycle'] = cycler(color=['blue','crimson','goldenrod','0.75'])

        cax = axs[v]

        #############################################################
        #### 24 hour cycle
        #############################################################
        if v == 0 :
            texts = 'a) 7-30nm'
        elif v == 1 :
            texts = 'b) 30-140nm'
        elif v == 2 :
            texts = 'c) 140-800nm'
        elif v == 3 :
            texts = 'd) 800nm+'                
 
       
        find_freq = np.abs(pr_dys-1.0).argmin() + 1
        print(str(pr_dys[find_freq-1])+' day cycle')
        hours_for_analysis = pr_dys[find_freq-1]*24;
        
        #pwr_thresh = 0.05
        
        x_utc = np.arange(0,25,2) # local daylight time
        x_utc = np.insert(x_utc,13,24.1)
        x_cdt = np.arange(0,25,2)-5 # utc time
        for i in np.arange(0,len(x_cdt)):
            if x_cdt[i] < 0:
                x_cdt[i] = x_cdt[i]+24

        plt.rcParams.update({'font.size': 14});

        cnt_plt = 0; max_val = 0
        for i in np.arange(1,5):
            
            # Choose power thresholdhold to be greater than the 99th percentile above red noise for each individual Season
            pwr_thresh = copy.deepcopy(F99red2[i][find_freq-1])
            pwr_thresh = copy.deepcopy(red_noise_norm[i,find_freq-1])

            print(pwr_thresh)
            
            temp_times = chunktime_s_comb[i]
            temp_data = copy.deepcopy(maxt_s_comb[i][find_freq,:])
            temp_pwr_data = copy.deepcopy(pwr_s_comb[i][:,find_freq])
            for ii in np.arange(0,len(temp_data)):
                if temp_pwr_data[ii] < pwr_thresh or np.isnan(temp_pwr_data[ii]):
                    temp_data[ii] = np.nan
                    
            cnt_vals = np.sum(np.logical_not(np.isnan(temp_data)))
            temp_data = temp_data[np.logical_not(np.isnan(temp_data))]
        
            bins = np.arange(0,25,2) 
            a1,a2 = np.histogram(temp_data,bins=bins)
            a1 = np.insert(a1,0,0)
            a1 = np.insert(a1,13,0)

            a1_n = a1/np.sum(a1)

            
        
            cax.plot(x_utc,a1_n,ls='steps',lw=5-(0.5*i),color=colors[i],label=cnames_plt[cnt_plt]+' ('+str(cnt_vals)+')')
            cax.set_xticks(x_utc)
        
            cnt_plt = cnt_plt + 1
        
            max_val = np.max([np.max(a1_n),max_val])
        
        #cax.set_ylim([0,max_val+0.05])        
        #cax.text(0,max_val-0.05,texts,fontsize = 16)

        cax.set_ylim([0,0.7])        
        cax.text(0,0.6,texts,fontsize = 15)

        cax.set_ylabel('Normalized Frequency')

        if v == 0 or v == 3:
            cax.legend(fontsize=13,loc="upper center")
        else:
            cax.legend(fontsize=13,loc="upper right")

        if v < 3:
            cax.tick_params(axis='x',labelbottom='off')

        if v == 3 :
            cax.set_xlabel('Hour of Day (UTC)')

            plt.rcParams['xtick.labelsize']=14
            ax2 = cax.twiny()
            ax2.set_xticks(x_utc)
            ax2.set_xticklabels(x_cdt,fontsize=14)
            ax2.xaxis.set_ticks_position('bottom')
            ax2.xaxis.set_label_position('bottom')
            ax2.spines['bottom'].set_position(('outward',50))
            ax2.set_xlabel('Hour of Day (Central Daylight Time)',fontsize=14)
            ax2.set_xlim(cax.get_xlim())

    #plt.legend(fontsize=12)
if plot_combined == 1:
    plt.tight_layout()
    plt.savefig(savepath+addsave+'Power_Spectrum_1FIG_Separate.png')
    plt.savefig(savepath+addsave+'Power_Spectrum_1FIG_Separate.eps')
    #plt.close(fig)
if plot_combined == 2:
    plt.tight_layout()
    plt.savefig(savepath+addsave+'24Hour_Cycle_Timing_1FIG.png')
    plt.savefig(savepath+addsave+'24Hour_Cycle_Timing_1FIG.eps')
    plt.close(fig)    