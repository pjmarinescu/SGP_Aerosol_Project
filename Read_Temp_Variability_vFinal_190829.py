#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:09:58 2017

@author: pmarin
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
import pandas as pd
import copy
import glob
import scipy

######################################
### READ in raw data from Ezra Files
######################################
final_path = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/RP_072720_i/'
#rootdir = final_path+'Data/processed/'
#rootdir = '/Volumes/Data/Act_Spectrum/processed/'
#rootdir = '/Users/pmarin/Desktop/processed/' # Downloaded from above /Volume/Data/Act_Spectrum

#rootdir = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/processed/'
rootdir = '/Users/pmarin/Documents/processed/' # Downloaded from above

nb = 1000; # number of bins in CCN_spectra file
nd = 215; # number of bins in size distribution

# Read in Activatation Spectrum Data From Files fro Ezra Levin
act_data = np.zeros((1,nb)); intN = np.zeros((1,1)); cnt = 0; 
dist = np.zeros((1,nd)); ave_k = np.zeros((1,1))
distDP = np.zeros((1,nd)); kdist = np.zeros((1,215))
#ccn_data = np.zeros((1,5)); 
ccn_data = [0]; ccn_ss = [0]; cpc_data = np.zeros((1,1))
cpc_cnt =[];times = [0]

for subdir, dir, files in os.walk(rootdir):
    cur_dir = os.path.join(subdir); print(cur_dir)

    f = open('/Volumes/Data/Act_Spectrum/SGP_CCN/Feb2018_FromEzra/20110101.080840/Dp.txt', 'r')
    standDP = np.array(f.read().split(), dtype=float)
    f.close()
    
    fileDP = os.path.join(subdir,"Dp.txt")    
    fileTime = os.path.join(subdir,"time.txt")
    fileSpec = os.path.join(subdir,"CCN_spectra.txt")
    fileIntn = os.path.join(subdir,"integrated_N_tot.txt")
    fileDist = os.path.join(subdir,"dndlogdp.txt")
    filek = os.path.join(subdir,"ave_k.txt")
    filekDist = os.path.join(subdir,"k_dist.txt")
    fileDim = os.path.join(subdir,"dimensions.txt")
    fileCPC = os.path.join(subdir,"Measured_N_conc.txt")
    fileCCN = os.path.join(subdir,"Measured_ccn.txt")
    fileCCNset = os.path.join(subdir,"Measured_ccn_ssw.txt")

    # file counter (ignore .DS.tore)
    if cnt == 0:
        cnt = cnt + 1
        continue

    # Read time from tile
    with open(fileTime) as f:
        time_temp = np.array(f.read().split(), dtype=float)   
    nt = time_temp.size # Get time dimension from number of points in time file

    if time_temp[nt-1] == 0:
        time_temp[nt-1] = 86400

    # Skipe Data when times are not in ascending order -- often if two instances with the same times
    if any(n < 0 for n in np.diff(time_temp)):
        print('Bad Times for:'+fileTime)
        continue

    # Read time from filenames
    time_string = subdir[len(subdir)-8:len(subdir)]
    file_date = datetime.strptime(time_string,'%Y%m%d')
    print(file_date)

    # Get dimensions of files
    f = open(fileDim, 'r')
    dims = np.array(f.read().split(), dtype=float)
    f.close()

    # Get diameter dimension for size distribution
    f = open(fileDP, 'r')
    DP = np.array(f.read().split(), dtype=float)
    f.close()

    DP2 = np.zeros((1,nd)); DP2[:] = DP
    
    for t in np.arange(0,nt):
        distDP = np.append(distDP,DP2,axis=0)

    # Read particle size distributions
    with open(fileDist) as f:
        dist_temp = np.array(f.read().split(), dtype=float).reshape(nt, int(dims[3]))
        #DP_cur = standDP
        dist = np.append(dist,dist_temp,axis=0)

        del dist_temp
        f.close()       
    
    datetime_object = datetime.strptime(time_string[0:8], '%Y%m%d') 
    for i in np.arange(0,nt,1):
        times.append(datetime_object + timedelta(0,time_temp[i]))     

    # Read Calculated CCN Spectra   
#    with open(fileSpec) as f:
#        #print(np.array(f.read().split(), dtype=float).reshape(nt, nb))
#        print(np.shape(np.array(f.read().split(), dtype=float)))
#        data_temp = np.array(f.read().split(), dtype=float).reshape(nt, nb)
#        act_data = np.append(act_data,data_temp,axis=0)
#        del data_temp
#        f.close()

    # Read Integrated Number from SMPS/DMA  
    with open(fileIntn) as f:
        intN_temp = np.array(f.read().split(), dtype=float).reshape(nt, 1)
        intN = np.append(intN,intN_temp,axis=0)
        del intN_temp
        f.close()

    # Read average kappa value from tdma
    with open(filek) as f:
        k_temp = np.array(f.read().split(), dtype=float).reshape(nt, 1)
        ave_k = np.append(ave_k,k_temp,axis=0)
        del k_temp
        f.close()

    # Read kappa distribution from tdma
    #with open(filekDist) as f:
    #    kdist_temp = np.array(f.read().split(), dtype=float).reshape(nt, int(dims[3]))
    #    kdist = np.append(kdist,kdist_temp,axis=0)
    #    del kdist_temp
    #    f.close()

    # Read integrated number from cpc  
    if os.path.isfile(fileCPC):
        cpc_cnt.append(cnt)
        with open(fileCPC) as f:
            cpc_temp = np.array(f.read().split(), dtype=float).reshape(nt, 1)
            cpc_data = np.append(cpc_data,cpc_temp,axis=0)
            del cpc_temp
            f.close()
    else:
        cpc_cnt.append(cnt)
        cpc_temp = np.zeros((nt,1))
        cpc_temp[:] = np.nan
        cpc_data = np.append(cpc_data,cpc_temp,axis=0)
        del cpc_temp

    # Read ccn spectra from ccn counter 
    if os.path.isfile(fileCCN) and dims[4] != 1:
        with open(fileCCN) as f:
            ccn_temp = np.array(f.read().split(), dtype=float).reshape(nt, int(dims[4]))
            for t in np.arange(0,nt):
                ccn_data.append(ccn_temp[t,:])

            del ccn_temp
            f.close()
    else:
        ccn_temp = np.zeros((nt,7))
        ccn_temp[:] = np.nan
        for t in np.arange(0,nt):
            ccn_data.append(ccn_temp[t,:])
        del ccn_temp
            
    # Read ccn set supersaturations from ccn counter 
    if os.path.isfile(fileCCNset) and dims[4] != 1:
        with open(fileCCNset) as f:
            ccn_temp = np.array(f.read().split(), dtype=float).reshape(nt, int(dims[4]))
            for t in np.arange(0,nt):
                ccn_ss.append(ccn_temp[t,:])
            del ccn_temp
            f.close()
    else:
        ccn_temp = np.zeros((nt,7))
        ccn_temp[:] = np.nan
        for t in np.arange(0,nt):
            ccn_ss.append(ccn_temp[t,:])
        del ccn_temp   
    cnt = cnt + 1


# convert bin midpoints to bin edges and dlogDP
standDP_edges = (standDP[0:len(standDP)-1]+standDP[1:])/2
standDP_edges = np.concatenate((standDP_edges[0]-(np.diff(standDP_edges[0:2])),standDP_edges,standDP_edges[-1]+(np.diff(standDP_edges[-2:]))),axis=0)
dlogDP = np.diff(np.log(standDP_edges))
dlog10DP = np.diff(np.log10(standDP_edges))
dist = dist*np.log10(np.exp(1)) # Convert everything to natural log


# Load in Peter's Corrected CPC data
pd_cpc = pd.read_excel(final_path+'CPC5.xlsx')
cpc_data5 = pd_cpc.iloc[:,0].values
cpc_data5 = np.insert(cpc_data5,0,0)

#####################################################################
#############
#############   Exclude Bad Data Points (QCing)
#############
#####################################################################
# Create QCing masks that can be applied to the data
n_masks = 8;
mask_cpcltdma = np.zeros(len(times))
mask_cpc0 = np.zeros(len(times))
mask_both = np.zeros(len(times))
mask_intNbad = np.zeros(len(times)) 
mask_negk = np.zeros(len(times))
mask_noact = np.zeros(len(times))
mask_nanpsd = np.zeros(len(times))
mask_psd = np.zeros(len(times))
mask_cal = np.zeros(len(times))
mask_data = np.zeros((len(times),n_masks))

##################################################
# Define Variables 
intNcalc = np.zeros(np.shape(intN))
intdiff = np.zeros(np.shape(intN))
gact = np.zeros(np.shape(act_data));
gcpc = np.zeros(np.shape(cpc_data));
gafrac = np.zeros(np.shape(act_data)); 
gdist = np.zeros(np.shape(dist))
gkdist = np.zeros(np.shape(kdist))
gdist2 = np.zeros(np.shape(dist))
gdfrac = np.zeros(np.shape(dist)); gintN = np.zeros(np.shape(intN)) ; 
gDP = np.zeros(np.shape(distDP));

gccn = []; gccn_ss = []
bccn = []; bccn_ss = []
bact = np.zeros(np.shape(act_data));

b_dist_nanpsd = np.zeros(np.shape(dist)); 
b_dist_psd = np.zeros(np.shape(dist))
b_intN = np.zeros(len(times))
b_k = np.zeros(len(times))
b_kdist = np.zeros(np.shape(kdist))
b_noact = np.zeros(np.shape(act_data))

cpcLT = np.zeros((len(times),2))
cnt1 = 0; cntint = 0; cntnegk = 0; cntcpc0 = 0; cntcpclt = 0; cntnoact = 0; cntnanpsd = 0; cntpsd = 0; cntcal = 0;
gtimes = []; btimes0 = []; btimes1 = []; btimes2 = []; btimes3 = []; btimes4 = []; btimes5 = []; btimes6 = []; btimes7 = [];
cnt_dif1 = 0
dist_0 = copy.deepcopy(dist);
dist_0[np.isnan(dist_0)] = 0 
##################################################

# Note count begins at 1 [because read in arrays above were initialized with a zero)
for i in np.arange(1,len(times)):
    intNcalc[i] = np.sum(dist_0[i,:]*dlogDP) # 
    intdiff[i] = intNcalc[i] - intN[i]

    # Bad cpc data
    if cpc_data5[i] < 100 or np.isnan(cpc_data5[i]):
        #print("CPC = 0,NaN",times[i].strftime('%Y%m%d.%H%M%S'))  
        mask_cpc0[i] = 1
        mask_data[i,0] = 1
        btimes0.append(times[i])
        cntcpc0 = cntcpc0 + 1

    # Bad Integrated Number (unrealistically high or low)
    if np.isnan(intN[i]) or (intN[i] > 100000) or (intN[i] < 100): #intN[i] == 0 or
        #print("Bad Integrated Number at cnt:",times[i].strftime('%Y%m%d.%H%M%S'))  
        b_intN[cntint] = intN[i]
        btimes2.append(times[i])
        mask_intNbad[i] = 1
        mask_data[i,2] = 1
        cntint = cntint + 1

#    # Bad Kappa Values in Kappa Distributions
#    if np.sum(n < 0 for n in kdist[i,:]) > 0:
#        #print("Negative Kappa:",times[i].strftime('%Y%m%d.%H%M%S'))
#        b_k[cntnegk] = ave_k[i]        
#        b_kdist[cntnegk,:] = kdist[i,:]        
#        mask_negk[i] = 1
#        #mask_data[i,3] = 1
#        btimes3.append(times[i])
#        cntnegk = cntnegk + 1

#    # No activation spectrum data
#    if np.nansum(act_data[i,:]) == 0:
#        #print("No Activation Spectrum Data at cnt:",times[i].strftime('%Y%m%d.%H%M%S'))
#        mask_noact[i] = 1
#        mask_data[i,4] = 1
#        btimes4.append(times[i])
#        b_noact[cntnoact,:] = act_data[i,:] 
#        cntnoact = cntnoact + 1


    # Missing supermicron particles (i.e. APS data)       
#    if np.nansum(dist[i,127:]) == 0:
#        #print("NaN in PSD Data at small sizes (between 12 and 30 nm):",times[i].strftime('%Y%m%d.%H%M%S'))  
#        mask_noact[i] = 1
#        mask_data[i,4] = 1
#        btimes4.append(times[i])
#        cntnoact = cntnoact + 1

        
    # Distribution shifted, so no longer capturing smallest particles       
    if np.nansum(dist[i,0:23]) == 0:
        #print("NaN in PSD Data at small sizes (between 12 and 30 nm):",times[i].strftime('%Y%m%d.%H%M%S'))  
        mask_nanpsd[i] = 1
        mask_data[i,5] = 1
        btimes5.append(times[i])
        b_dist_nanpsd[cntnanpsd,:] = dist[i,:]
        cntnanpsd = cntnanpsd + 1

    # Unrealistic Peaks in the size Distribution        
    if (np.nanmax(dist[i,15:] / intN[i]) > 3) or (np.nanmax(dist[i,100:] / intN[i]) > 0.5) :
        #print('Unrealastic Peaks in Size Dist:',times[i].strftime('%Y%m%d.%H%M%S'))
        mask_psd[i] = 1
        mask_data[i,6] = 1
        btimes6.append(times[i])
        b_dist_psd[cntpsd,:] = dist[i,:] / intN[i]
        cntpsd = cntpsd + 1
    
    # Wierd peaks right in intN right after DMA calibrated (i.e. after a long time break)        
    if (times[i].hour < 10) and (times[i].hour > 7):
        cpc_pct_b = (cpc_data5[i]-cpc_data5[i-1])/cpc_data5[i-1]
        int_pct_b = (intN[i,:]-intN[i-1,:])/intN[i-1,:]        
        cpc_pct_f = (cpc_data5[i+1]-cpc_data5[i])/cpc_data5[i]
        int_pct_f = (intN[i+1,:]-intN[i,:])/intN[i,:]        
#        if (int_pct > 3 * cpc_pct) and (int_pct - cpc_pct) > 0.25:
        if (int_pct_b - cpc_pct_b) > 0.10 and (cpc_pct_f - int_pct_f) > 0.10 and (times[i]-times[i-1]).seconds/60 > 100: #100 minutes
            mask_cal[i] = 1
            mask_data[i,7] = 1
            btimes7.append(times[i])
            cntcal = cntcal + 1    
        elif (int_pct_b - cpc_pct_b) > 0.50 and (cpc_pct_f - int_pct_f) > 0.50:
            mask_cal[i] = 1
            mask_data[i,7] = 1
            btimes7.append(times[i])
            cntcal = cntcal + 1 
    
    # If passes tests, data is good (g)
    if np.sum(mask_data[i,:]) == 0:

        gccn.append(ccn_data[i]) 
        gccn_ss.append(ccn_ss[i]) 
        gcpc[cnt1,:] = copy.deepcopy(cpc_data5[i])
        #gact[cnt1,:] = copy.deepcopy(act_data[i,:]) 
#        gafrac[cnt1,:] = act_data[i,:] / intN[i]
        #gkdist[cnt1,:] = copy.deepcopy(kdist[i,:])
        gdist[cnt1,:] = copy.deepcopy(dist[i,:])
        gdfrac[cnt1,:] = dist[i,:] / intN[i]
        gDP[cnt1,:] = copy.deepcopy(distDP[i,:])
        gintN[cnt1,:] = copy.deepcopy(intN[i,:])
        gtimes.append(times[i])

#        print(np.nansum(gdist2[cnt1,:-1]*np.diff(np.log10(gDP[cnt1,:]))),cpc_data[i,:],intN[i,:])   
#        if np.abs((np.nansum(gdist2[cnt1,:-1]*np.diff(np.log10(gDP[cnt1,:]))) - cpc_data[i,:])/cpc_data[i,:]) > 0.1:
#            break    
            
        cnt1 = cnt1 + 1

mask_data_sum = np.sum(mask_data,axis=1)
total_eliminated_points = np.sum(mask_data_sum>0); 
print('Good:'+str(cnt1)+' | Eliminated:'+str(total_eliminated_points))
gtimes_arr = np.asarray(gtimes)


cpc_bad_flag = np.zeros(cnt1+1)
for i in np.arange(0,cnt1):
    if gcpc[i,:] < 100 or np.isnan(gcpc[i,:]):
        cpc_bad_flag[i] = 1

######## Process met dataframe
def standardize(x):
    y = (x-np.nanmean(x))/np.nanstd(x)
    return y

###############################
###### Plot Qcing Flags with Concentrations
################################
#import matplotlib   
#import matplotlib.dates as mdates
#for yyyy in np.arange(2009,2014,1):
#    for t in np.arange(1,13,1):
#        
#        if t == 2:
#            dd = 28
#        else:
#            dd = 30
#        
#        xlims = matplotlib.dates.date2num([datetime(yyyy,t,1), datetime(yyyy,t,dd)])
#        
#        fig,ax = plt.subplots(1,1,figsize=(10,5))
#        plt.plot(np.array(times[1:]),mask_data[1:,0]*1,'ok',label='CPC')
#        plt.plot(np.array(times[1:]),mask_data[1:,2]*2,'^k',label='intN')
#        plt.plot(np.array(times[1:]),mask_data[1:,5]*3,'dk',label='dist')
#        plt.plot(np.array(times[1:]),mask_data[1:,6]*4,'xk',label='peaks')
#        plt.plot(np.array(times[1:]),mask_data[1:,7]*5,'*k',label='Cal')    
#        ax.set_ylabel('QC Flags')
#        plt.title(str(t)+' / '+str(yyyy))
#        plt.legend(fontsize=10,loc='upper left')
#        
#        if t < 10:
#            strt = '0'+str(t)
#        else:
#            strt = str(t)
#        
#        ax2 = ax.twinx()
#        ax2.plot(df_cpc.index,df_cpc.cpc_avg,'og',ms=2,label='Full CPC Data')
#        ax2.plot(np.array(times[1:]),cpc_data[1:],'ob',ms=1.5,label='Ezra CPC Data')
#        ax2.plot(np.array(times[1:]),intN[1:],'om',ms=1,label='IntN Data')
#        ax2.plot(np.array(times[1:]),cpc_data5[1:],'oy',ms=1,label='Peter CPC Data')
#        ax2.set_xlim(xlims)
#        ax2.set_ylim([50,50000])
#        ax2.set_ylabel('CPC Conc. (#/cc)')
#        ax2.set_yscale('log')
#        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
#        plt.legend(fontsize=10,loc='upper right')
#        plt.tight_layout()
#        plt.savefig('/Volumes/Data/Act_Spectrum/QC_Temp_Var_'+str(yyyy)+strt+'_v0828.png')
#        plt.close(fig)
#        


#####################################################################
#############
#############   Load CCN Data from Ezra (QCing)
#############   Take into account CCNC data for adjusting size distributions
#############   Not include in V1 of size distribution product, including in V2
#####################################################################

pathname_ccn = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/ccn_ratio_RP072120_noi/'
pathname_ccn = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/ccn_ratio_RP072120_i/'
ccn_files = glob.glob(pathname_ccn+'*')

ccn_dates = []; ccn_ratio_SA = []; ccn_ratio_CSA = [] # Initialize arrays for data
for i in np.arange(0,len(ccn_files)):
    file_ccn = ccn_files[i]

    # Read data from Ezra Files    
    f = open(file_ccn, 'r')
    data_ccn_E = np.array(f.read().split(), dtype=float)
    data_ccn_E = np.reshape(data_ccn_E,[int(len(data_ccn_E)/3), 3])
    f.close()
    

#    yyyy = int(file_ccn[135:139])
#    mm = int(file_ccn[139:141])
#    dd = int(file_ccn[141:143])

#   RP_i
    yyyy = int(file_ccn[137:141])
    mm = int(file_ccn[141:143])
    dd = int(file_ccn[143:145])

#   RP_noi
#    yyyy = int(file_ccn[139:143])
#    mm = int(file_ccn[143:145])
#    dd = int(file_ccn[145:147])

    date_init = datetime(yyyy,mm,dd)
    
    for t in np.arange(0,np.shape(data_ccn_E)[0]):
        ccn_dates.append(date_init+timedelta(seconds=int(data_ccn_E[t,0])))
        ccn_ratio_SA.append(data_ccn_E[t,1])
        ccn_ratio_CSA.append(data_ccn_E[t,2])        

ccn_dates_arr = np.asarray(ccn_dates)
ccn_ratio_SA = np.asarray(ccn_ratio_SA)  # CCN spectra / SA estimates
ccn_ratio_CSA = np.asarray(ccn_ratio_CSA) # CCN spectra / CSA estimates
ccn_ratio_SA[ccn_ratio_SA == 0.0] = np.nan 
ccn_ratio_CSA[ccn_ratio_CSA == 0.0] = np.nan 

# Bad CCN data: https://www.archive.arm.gov/ArchiveServices/DQRService?dqrid=D141208.1
ccn_id_bad1 = np.abs(ccn_dates_arr-datetime(2009,5,15,12)).argmin() ; 
ccn_id_bad2 = np.abs(ccn_dates_arr-datetime(2009,8,3,0)).argmin()
ccn_ratio_SA[ccn_id_bad1:ccn_id_bad2+1] = np.nan
ccn_ratio_CSA[ccn_id_bad1:ccn_id_bad2+1] = np.nan

# shows too low compared to SMPS, while NPF compaign showed no biases
ccn_id_bad1 = np.abs(ccn_dates_arr-datetime(2013,2,27,0)).argmin() ; 
ccn_id_bad2 = np.abs(ccn_dates_arr-datetime(2013,6,15,0)).argmin()
ccn_ratio_SA[ccn_id_bad1:ccn_id_bad2+1] = np.nan
ccn_ratio_CSA[ccn_id_bad1:ccn_id_bad2+1] = np.nan

gtimes_sec = np.zeros(cnt1)
for i in np.arange(0,cnt1):
    gtimes_sec[i] = (gtimes_arr[i]-datetime(2009,1,1,0)).total_seconds()

ccn_sec = np.zeros(len(ccn_ratio_SA))
for i in np.arange(0,len(ccn_ratio_SA)):
    ccn_sec[i] = (ccn_dates_arr[i]-datetime(2009,1,1,0)).total_seconds()

#################### -------------------- FINISH READING IN DATA


###############################################################################
###############################################################################
######              Adjust Size Distributions of "Good" Data
###############################################################################
###############################################################################
gdist2 = copy.deepcopy(gdist)
# Corrections for Remaining Particles after systematic bias correction
pathname = copy.deepcopy(final_path)
savename = 'vFinal_D2E4_var99_LT2wk18'     # Extrapolation: Diameter^2 function
savename = 'vFinal_D2E4_var95_LT2wk18_CPCQC'               # Expontential Correction: exp(4) function
savename = 'vFinal_D2E4_var95_LT2wk18_CPCQC_medpctdif'     # Expontential Correction: exp(4) function
                                            # var99 -- only do exponential correction up to the 99th percentile diameter
                                            # Long term correction: 2wk median bias, only including hours 0-18

savename = 'vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072720_i' 

ccn_corr_flag = 1 # 0 if exclude, V1 dataset
                  # 1 if include, V2 dataset

#########################################################
#########################################################
# Extrapolate DMA size distribution down to these values
#########################################################
#########################################################
EXT_DP = copy.deepcopy(standDP)
logDP_stand = np.log(EXT_DP)
diam_low = 0.007 # smallest diameter to exrapolate down to

DP_stand_edges = (standDP[0:len(standDP)-1]+standDP[1:])/2
DP_stand_edges = np.concatenate((DP_stand_edges[0]-(np.diff(DP_stand_edges[0:2])),DP_stand_edges,DP_stand_edges[-1]+(np.diff(DP_stand_edges[-2:]))),axis=0)
logDP_edges = np.log(DP_stand_edges)
dlogDP_mids = np.diff(logDP_edges)
dlogDP_median = np.nanmedian(dlogDP_mids[0:10])

# Extend diameter size bins down to diam_low
logDP_extend = []; DP_extend = []
for i in np.arange(1,1000):
    logDP_extend.append(np.log(standDP[0])-(dlogDP_median*i))
    DP_extend.append(np.power(np.exp(1),logDP_extend[i-1]))
    if np.power(np.exp(1),logDP_extend[i-1]) < diam_low:
        break

# Caclculate diameter bins, midpoints, and differences
DP_extend = DP_extend[::-1]
DP_full = np.concatenate((DP_extend,EXT_DP),axis=0)
DP_full_mids = (DP_full[0:len(DP_full)-1]+DP_full[1:])/2
DP_full_mids = np.concatenate((DP_full_mids[0]-(np.diff(DP_full_mids[0:2])),DP_full_mids,DP_full_mids[-1]+(np.diff(DP_full_mids[-2:]))),axis=0)
dlogDP_full = np.diff(np.log(DP_full_mids))
logDP_full = np.log(DP_full)

# Actual extrapolation step
gdist_extrap = np.zeros((cnt1,len(DP_full)))
for i in np.arange(0,cnt1):         
    if cpc_bad_flag[i] == 1: # only extrapolate functions were CPC data available
        gdist_extrap[i,:] = np.nan    
    else:
        nonnan_id =  next(x for x in np.arange(0,len(gdist[i,:])) if not np.isnan(gdist[i,x])) # find first non nan value                                        
        # Define function to fit to data and for extrapolation
        def func(x,a,b):
            return a * np.power(x,2) + b # a * x^2 + b
        
    #    Tested additonal functions below, but the above function seemed to work best
    #    def func2(x,a,b,c):
    #        return a * np.power(x,2) + b*x + c # a * x^2 + b
    #    def func3(x,a):
    #        return a * np.power(b,np.log10(x)) # a * x^2 + b
    
        d_nbins = len(DP_full)-len(standDP) # number of additional bins due to the extrapolation
        z1, z2 = scipy.optimize.curve_fit(func,standDP[nonnan_id:nonnan_id+5],gdist[i,nonnan_id:nonnan_id+5]) #extrapolate using the data up to the 5th value after the first non-nan value    
        gdist_extrap[i,nonnan_id+d_nbins:] = copy.deepcopy(gdist[i,nonnan_id:]) # Replace all non-nan ideas with the original data
        gdist_extrap[i,:nonnan_id+d_nbins] = func(DP_full[:nonnan_id+d_nbins],*z1)
        gdist_extrap[gdist_extrap<0] = 0 # No negative values in size distribution    
        print(i)
#########################################################
gdist_extrap_save = copy.deepcopy(gdist_extrap) # Save values after Step 1: Extrapolation
#########################################################
#########################################################

gdist_extrap_cpc = copy.deepcopy(gdist_extrap)
gintN_full = np.zeros(cnt1)
cur_cpcint_dif = np.zeros(cnt1); 
cur_cpcint_pdif = np.zeros(cnt1);

#########################################################
#########################################################
# Adjust Extrapolated Distributions for "CPC Vision"
# Values from Mertes et al. 
#########################################################
#########################################################
DP_CPC = [0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.020,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.030]
EFF_CPC =[0.000,0.010,0.080,0.195,0.325,0.450,0.550,0.650,0.730,0.800,0.850,0.890,0.920,0.940,0.955,0.965,0.975,0.980,0.985,0.990,0.995,0.996,0.997,0.998,0.999,1.000]
for i in np.arange(0,cnt1,1):   
#for i in np.arange(idx_npf1,idx_npf2):   
    # Adjust cpc amount distributions for CPC efficiency       
    if cpc_bad_flag[i] == 1: # if no cpc data, don't do any of the processing
        gdist_extrap_cpc[i,:] = np.nan
        gintN_full[i] = np.nan
    else:        
        for k in np.arange(0,len(DP_full)-2):                
            if DP_full[k] > 0.031: #or k >= (nonnan_id): 
                break                    
            cpc_deteff = np.interp(DP_full[k],DP_CPC,EFF_CPC) # Detection efficient interpolation based on Mertes paper and dT = 17C
            gdist_extrap_cpc[i,k] = gdist_extrap_cpc[i,k] * cpc_deteff                  
        gintN_full[i] = np.nansum(gdist_extrap_cpc[i,:]*dlogDP_full)
########################################################
gdist_extrap_cpc_save = copy.deepcopy(gdist_extrap_cpc) # Save values after Step 2: Scaling Down for CPC Efficiencies

#########################################################
#########################################################
# Correct for Long-term bias (2-week)
#########################################################
#########################################################

# Calculate difference between CPC in new integrated DMA values post Steps 1 and 2
ltbias_gintN_full = np.zeros(cnt1)
ltbias_cpcint_dif = np.zeros(cnt1)
ltbias_cpcint_pdif = np.zeros(cnt1)
ltbias_cpcint_ratio = np.zeros(cnt1)
opp_ltbias_cpcint_pdif = np.zeros(cnt1)
cnt_no_reconcile = 0
for i in np.arange(0,cnt1):

    if cpc_bad_flag[i] == 1: # if no cpc data, don't do any of the processing
        ltbias_cpcint_dif[i] = np.nan
        ltbias_cpcint_pdif[i] = np.nan  
        ltbias_cpcint_ratio[i] = np.nan      
        opp_ltbias_cpcint_pdif[i] = np.nan
    else:
        ltbias_cpcint_dif[i] = gcpc[i] - gintN_full[i]
        ltbias_cpcint_pdif[i] = ltbias_cpcint_dif[i] / gintN_full[i]      
        ltbias_cpcint_ratio[i] = gcpc[i] / gintN_full[i]      
        opp_ltbias_cpcint_pdif[i] = (gintN_full[i] - gcpc[i]) / gcpc[i]      
        # Only include hours 0 through 18, since we would expect the two to be potentially different new particle formation times.
        if gtimes[i].hour >=18: 
            #gdist_extrap_cpc[i,:] = np.nan
            ltbias_cpcint_pdif[i] = np.nan
            ltbias_cpcint_dif[i] = np.nan
            ltbias_cpcint_ratio[i] = np.nan      
            cnt_no_reconcile = cnt_no_reconcile + 1 

ltbias_cpcint_pdif_meananom = ltbias_cpcint_pdif - np.nanmean(ltbias_cpcint_pdif)
ltbias_cpcint_pdif_medanom = ltbias_cpcint_pdif - np.nanmedian(ltbias_cpcint_pdif)

print(cnt_no_reconcile)

gdist_extrap_cpc_lt = np.zeros(np.shape(gdist_extrap_cpc))

lt_gcpc = np.zeros(cnt1)
avg_pct_dif = np.zeros(cnt1)

med_pct_dif = np.zeros(cnt1)
med_pct_dif_meananom = np.zeros(cnt1)
med_pct_dif_medanom = np.zeros(cnt1)

med_ratio = np.zeros(cnt1)
med_ratio_ccn = np.zeros(cnt1)
med_ratio_both = np.zeros(cnt1)
med_pct_dif_both = np.zeros(cnt1)

lt_scale_factor = np.zeros(cnt1)
gintN_full_2 = np.zeros(cnt1)
cur_cpcint_dif = np.zeros(cnt1)
cur_cpcint_pdif = np.zeros(cnt1)

#med_pct_dif_1M = np.zeros(cnt1); med_pct_dif_1W = np.zeros(cnt1); med_pct_dif_3W = np.zeros(cnt1); med_pct_dif_3M = np.zeros(cnt1)
for i in np.arange(0,cnt1):

    # First, Adjust entire distribution by a long term bias value (weekly)
    day_m1 = gtimes[i]-timedelta(hours=168); m1_id = np.abs(gtimes_arr-day_m1).argmin()
    day_p1 = gtimes[i]+timedelta(hours=168); p1_id = np.abs(gtimes_arr-day_p1).argmin()

    m2_id = np.abs(ccn_dates_arr-day_m1).argmin()
    p2_id = np.abs(ccn_dates_arr-day_p1).argmin()

    lt_gcpc[i] = copy.deepcopy(gcpc[i])

    if cpc_bad_flag[i] == 1: # if no cpc data, don't do any of the processing
        avg_pct_dif[i] = np.nan; med_pct_dif[i] = np.nan
        med_pct_dif_meananom[i] = np.nan; med_pct_dif_medanom[i] = np.nan
        lt_scale_factor[i] = np.nan; gdist_extrap_cpc_lt[i,:] = np.nan
        gintN_full_2[i] = np.nan; cur_cpcint_dif[i] = np.nan; cur_cpcint_pdif[i] = np.nan; 
        med_ratio[i] = np.nan;  med_ratio_ccn[i] = np.nan
    else:
        avg_pct_dif[i] = np.nanmean(ltbias_cpcint_pdif[m1_id:p1_id+1])
        med_pct_dif[i] = np.nanmedian(ltbias_cpcint_pdif[m1_id:p1_id+1])
        med_ratio[i] = np.nanmedian(ltbias_cpcint_ratio[m1_id:p1_id+1])
        med_pct_dif_meananom[i] = np.nanmedian(ltbias_cpcint_pdif_meananom[m1_id:p1_id+1])
        med_pct_dif_medanom[i] = np.nanmedian(ltbias_cpcint_pdif_medanom[m1_id:p1_id+1])

        if np.count_nonzero(~np.isnan(ccn_ratio_SA[m2_id:p2_id+1])) > 40:
            med_ratio_ccn[i] = np.nanmedian(ccn_ratio_SA[m2_id:p2_id+1])
        else:
            med_ratio_ccn[i] = np.nan

        med_ratio_both[i] = np.nanmean([med_ratio[i],med_ratio_ccn[i]])
        med_pct_dif_both[i] = med_ratio_both[i] - 1
        
        # Caclulate long-term scale factor
        if ccn_corr_flag == 0:
            lt_scale_factor[i] = copy.deepcopy(med_pct_dif_medanom[i])
            lt_scale_factor[i] = copy.deepcopy(med_pct_dif[i])
        elif ccn_corr_flag == 1:
            lt_scale_factor[i] = copy.deepcopy(med_pct_dif_both[i])
            print('Yes')
            
        
        # Scaling distribuiton up / down depending on lt_scale_factor
        lt_gcpc[i] = copy.deepcopy(gcpc[i])
        gdist_extrap_cpc_lt[i,:] = gdist_extrap_cpc[i,:] + (gdist_extrap_cpc[i,:] * lt_scale_factor[i]) 
        gintN_full_2[i] = np.nansum(gdist_extrap_cpc_lt[i,:]*dlogDP_full)
    
        cur_cpcint_dif[i] = lt_gcpc[i] - gintN_full_2[i]
        cur_cpcint_pdif[i] = cur_cpcint_dif[i] / gintN_full_2[i]       


import pickle
pickle.dump([lt_scale_factor,gintN_full_2,cur_cpcint_dif,cur_cpcint_pdif,lt_gcpc],open('/Users/pmarin/Desktop/'+savename+".p","wb"))
#asdfasdf

########################################################
# Save new distributions including long-term bias corrections
gdist_extrap_cpc_lt_save = copy.deepcopy(gdist_extrap_cpc_lt) # Save size distribution after scaling for Systematic, Long-term Bias Correction

### PLOTTING SECTION TO CHECK THE LONG-TERM BIAS CORRECTIONS (for testing)
# 
# Calculation of scaling factor for the entire time period
#np.nanmedian(med_pct_dif)
#med_pct_dif_elim = copy.deepcopy(med_pct_dif)
#for i in np.arange(0,len(gtimes)):
#    if gtimes[i] > datetime(2009,1,1,0) and gtimes[i] < datetime(2009,3,20,0):
#        med_pct_dif_elim[i] = np.nan    
#    if gtimes[i] > datetime(2013,9,1,0) and gtimes[i] < datetime(2014,1,1,0):
#        med_pct_dif_elim[i] = np.nan    
#np.nanmedian(med_pct_dif_elim)
#
#fig = plt.figure(figsize=(10,8))
##plt.plot(gtimes,med_pct_dif_medanom,'-c',label='Med_Pct_Dif_MedAnom_E5X')
#plt.scatter(gtimes,ltbias_cpcint_ratio,color='m',label='CPC/TDMA Ratio')
#plt.scatter(ccn_dates,ccn_ratio_SA,color='b',label='CCN/SA Ratio')
#plt.plot(gtimes,med_ratio,'-y',label='Med CPC/TDMA Ratio')
#plt.plot(gtimes,med_ratio_ccn,'-r',label='Med CCN/SA Ratio')
#plt.plot(gtimes,med_ratio_both,'-c',label='Med Both Ratio')
#
##plt.scatter(ccn_dates,ccn_ratio_CSA,label='CCN/CSA Ratio')
#plt.ylim([0,5])
#plt.legend()
#
np.nanmedian(lt_scale_factor)
#np.nanmedian(med_pct_dif_medanom)
fig = plt.figure(figsize=(10,8))
#plt.plot(gtimes,med_pct_dif_medanom,'-c',label='Med_Pct_Dif_MedAnom_E5X')
plt.plot(gtimes,med_pct_dif,'-b',label='Med_Pct_Dif_E5X')
plt.scatter(gtimes,cur_cpcint_pdif,label='Med_Pct_Dif_E5X')
#plt.plot(gtimes,med_pct_dif_medanom_lt18,'-m',label='Med_Pct_Dif_MedAnom_LT18')
#plt.plot(gtimes,med_pct_dif_lt18,'-r',label='Med_Pct_Dif_LT18')
plt.plot(gtimes,np.ones(len(gtimes))*0.5,'-k')
plt.plot(gtimes,np.ones(len(gtimes))*-0.5,'-k')
plt.plot(gtimes,np.zeros(len(gtimes)),'-k')
plt.legend()
#
#fig = plt.figure(figsize=(10,8))
#plt.plot(gtimes,gcpc[0:len(gtimes)],'-c')
#plt.plot(gtimes,gintN_full[0:len(gtimes)],'-r')
#plt.legend()

#########################################################################################################################################################
##########################################################################################################################################################

###############################################################
###############################################################
# Correct Remaining Difference using an Exponential Function
###############################################################
###############################################################
xs = np.arange(np.log10(DP_full[0]),np.log10(DP_full[len(DP_full)-1]),0.001)
ys = np.zeros(len(xs))    

ys2 = np.zeros(len(xs))    
ys_do = np.zeros(len(xs))    

for l in np.arange(0,len(xs)):       
    ys[l] = 0.01 * np.exp(4*(-xs[l]+xs[0])) 
#    if np.power(10,xs[l]) > 0.280: # 0.280 is the approximate value of the 99th percentile diameter (280 nm)
    if np.power(10,xs[l]) > 0.200: # 0.200 is the approximate value of the 95th percentile diameter (280 nm)
        ys[l] = 0
        
ys_temp = copy.deepcopy(ys)
ys_temp[ys_temp==0] = np.nan
ys_subtract = np.nanmin(ys_temp)
ys_temp = ys_temp - ys_subtract
ys_temp[np.isnan(ys_temp)]=0
ys = copy.deepcopy(ys_temp)
        
DDP = np.diff(DP_full_mids)
DDP_scale = np.zeros(len(DP_full))
for l in np.arange(0,len(DP_full)):
    DDP_scale[l] = DDP[l] / np.nanmax(DDP)
    
    DP_full_ys = np.zeros(len(DP_full))
    DP_full_ys2 = np.zeros(len(DP_full))
    for l in np.arange(0,len(DP_full)):       
        DP_full_ys[l] = np.interp(np.log10(DP_full[l]),xs,ys)
        DP_full_ys2[l] = np.interp(np.log10(DP_full[l]),xs,ys2)

    cpc_deteff_k = np.zeros(len(DP_full))
    for l in np.arange(0,len(DP_full)):       
        cpc_deteff_k[l] = np.interp(DP_full[l],DP_CPC,EFF_CPC)

plt.rcParams.update({'font.size': 12})
fig,axs = plt.subplots(1,1,figsize=(7,4))
l1 = axs.plot(DP_full,DP_full_ys*cpc_deteff_k,'-k',lw=6,label='Final (Exponential * CPC Detection)')            
l2 = axs.plot(DP_full,DP_full_ys,':',color='gray',lw=4,label='Exponential Correction')            
axs.set_ylabel('Fraction of Remaining \n Particles to Add/Subtract')            
axs.set_xscale('log')            

ax2 = axs.twinx()  
l3 = ax2.plot(DP_full,cpc_deteff_k,':m',lw=4,label='CPC Detection Efficiency')            
ax2.set_xscale('log')            
ax2.set_ylabel('Fraction of Particles \n Detected by CPC')            
axs.set_xlabel('Diameter ($\mu$m)')
ax2.spines['right'].set_color('m')
ax2.yaxis.label.set_color('m')
ax2.tick_params(axis='y',colors='m')

ls = l1+l2+l3
labs = [l.get_label() for l in ls]
plt.legend(ls,labs,loc="center right") 

plt.tight_layout()
plt.savefig(pathname+'Correction_Curves_'+savename+'_CurveFitDP2.png')
plt.savefig(pathname+'Correction_Curves_'+savename+'_CurveFitDP2.eps')
plt.close(fig)
       
#############################################################################################################
###################### Determine remaining difference between instruments and adjust the smallest sizes
#############################################################################################################
gdist2_intN = np.zeros(cnt1); gdist2_intN_est = np.zeros(cnt1)
int_dif = np.zeros(cnt1); int_dif_est = np.zeros(cnt1)
gdist_extrap_cpc_lt_exp = np.zeros(np.shape(gdist_extrap_cpc))
gdist_extrap_cpc_est = np.zeros(np.shape(gdist_extrap_cpc))
cnt_iterate = 0
diam_90 = np.zeros(cnt1); diam_95 = np.zeros(cnt1); diam_99 = np.zeros(cnt1)
for i in np.arange(0,cnt1,1):  
    if cpc_bad_flag[i] == 1:    # No CPC data, so cannot do this processing
        diam_90[i] = np.nan
        diam_95[i] = np.nan
        diam_99[i] = np.nan
        gdist_extrap_cpc_lt_exp[i,:] = np.nan
        gdist2_intN[i] = np.nan
        gdist_extrap_cpc_est[i,:] = np.nan          
        gdist2_intN_est[i] = np.nan
        int_dif[i] = np.nan
        int_dif_est[i] = np.nan

    else:
        cur_dist = gdist_extrap_cpc_lt[i,:] * dlogDP_full   # Convert dN / dlogDP into bins    
        gdist_extrap_cpc_lt_exp[i,:] = copy.deepcopy(gdist_extrap_cpc_lt[i,:])
    
        for b in np.arange(0,len(cur_dist)):
            cur_bin_frac =  np.nansum(cur_dist[0:b]) / np.nansum(cur_dist)
            # Determine the diameter at which X% of the number concentrations falls below
            if (cur_bin_frac > 0.897) & (cur_bin_frac < 0.907):
                diam_90[i] = DP_full[b] 
            if (cur_bin_frac > 0.95) & (cur_bin_frac < 0.96):
                diam_95[i] = DP_full[b] 
            if cur_bin_frac > 0.99:
                diam_99[i] = DP_full[b] 
                break
            
        for l in np.arange(0,len(xs)):       
            ys[l] = 0.01 * np.exp(4*(-xs[l]+xs[0]))    #### NEED TO change ys here!!!!
            if np.power(10,xs[l]) > diam_95[i]: # Only correct data up to the 99th percentile diameter (avoids changing regions of size distribution with very small concentrations)
                ys[l] = 0
    
        # IF MORE PARTICLES in INTEGRATED DMA VS. CPC, REMOVE PARTICLES
        # (CPC - DMA)
        if cur_cpcint_dif[i] < 0 and np.abs(cur_cpcint_pdif[i]) > 0.01:                
    
            for kk in np.arange(0,1000): # iterations
                for k in np.arange(0,len(DP_full)): # up to 30 nm
                     cur_dist[k] = cur_dist[k] - (np.abs(cur_cpcint_dif[i]) * np.interp(np.log10(DP_full[k]),xs,ys) * np.interp(DP_full[k],DP_CPC,EFF_CPC))
                     cur_dist[cur_dist < 0] = 0            
                temp_cpcint_dif = lt_gcpc[i] - np.nansum(cur_dist)
            
                if np.abs(temp_cpcint_dif/cur_cpcint_dif[i]) < 0.01 or temp_cpcint_dif >= 0:
                    break
    
            if kk == 999:
                cur_dist[:] = np.nan            
                cnt_iterate = cnt_iterate + 1
                print('no solution')
    
            gdist_extrap_cpc_lt_exp[i,:] = cur_dist / dlogDP_full
            gdist2_intN[i] = np.nansum(gdist_extrap_cpc_lt_exp[i,:]*dlogDP_full)
                
        # IF LESS PARTICLES in INTEGRATED DMA VS. CPC, ADD PARTICLES
        if cur_cpcint_dif[i] >= 0 and np.abs(cur_cpcint_pdif[i]) > 0.01:    
            
            for kk in np.arange(0,1000):
                for k in np.arange(0,len(DP_full)):
                    cur_dist[k] = cur_dist[k] + (np.abs(cur_cpcint_dif[i]) * np.interp(np.log10(DP_full[k]),xs,ys) * np.interp(DP_full[k],DP_CPC,EFF_CPC))
    #            cur_dist[0:43] = smooth(cur_dist[0:43],window_len = 4)
                temp_cpcint_dif = lt_gcpc[i] - np.nansum(cur_dist)
                
                if np.abs(temp_cpcint_dif/cur_cpcint_dif[i]) < 0.01 or temp_cpcint_dif <= 0:
                    break
            if kk == 999:
                cur_dist[:] = np.nan            
                cnt_iterate = cnt_iterate + 1
                print('no solution')
    
            gdist_extrap_cpc_lt_exp[i,:] = cur_dist / dlogDP_full                       
            gdist2_intN[i] = np.nansum(gdist_extrap_cpc_lt_exp[i,:]*dlogDP_full)
        
        del(cur_dist)
        ######################################################################
        ################### SCALE UP SIZE DISTRIBUTIONS for CPC Vision
        ######################################################################
        gdist_extrap_cpc_est[i,:] = copy.deepcopy(gdist_extrap_cpc_lt_exp[i,:])
        for k in np.arange(0,100):                
            if DP_full[k] > 0.031:
                break                    
            
            cpc_deteff = np.interp(DP_full[k],DP_CPC,EFF_CPC) # Detection efficient interpolation based on Mertes paper and dT = 17C
            #print(cpc_deteff)
            # scale bin dN by the cpc efficiency and calculate the new integrated number for that bin
            gdist_extrap_cpc_est[i,k] = gdist_extrap_cpc_lt_exp[i,k] / cpc_deteff              
    
        gdist2_intN_est[i] = np.nansum(gdist_extrap_cpc_est[i,:]*dlogDP_full)
    
        int_dif[i] = gdist2_intN[i]-lt_gcpc[i]; #print('Integrated Dif:'+str(int_dif[i])+ '|'+str(cur_cpcint_dif[i]))
        int_dif_est[i] = gdist2_intN_est[i]-lt_gcpc[i]; #print('Integrated Dif:'+str(int_dif_est[i])+ '|'+str(cur_cpcint_dif[i]))
        print(i)
    
# Save final distributions 
gdist_extrap_cpc_last_save = copy.deepcopy(gdist_extrap_cpc_lt_exp) # pre final scale-up for CPC detection efficiencies
gdist_extrap_cpc_last_est_save = copy.deepcopy(gdist_extrap_cpc_est)  # post final scale-up for CPC detection efficiencies

#######################################
####### PLOT Sample distribution Demonstrating Corrections
#######################################
strs = ['(a) ','(b) ','(c) ']
example_times=[]

example_times.append(datetime(2013,5,11,1,0))
#example_times.append(datetime(2013,4,29,7,0))
example_times.append(datetime(2013,4,29,11,0))
example_times.append(datetime(2013,4,21,19,10))

#timewant = datetime(2012,2,2,1,0)
#ex_id = np.abs(np.array(gtimes[0:cnt1])-timewant).argmin()

#fig,ax = plt.subplots(1,1,figsize=(7,4));
#ax.plot(standDP[3:],gdist[ex_id,3:],'-k',lw=5)
#ax.set_xscale('log')
#ax.set_ylabel('dN (dln $D_p$)$^{-1}$')
#ax.set_xlim([0.007,14])
#ax.set_ylim([0,3000])
#ax.set_xlabel('Particle Diameter ($D_p$, $\mu$m)')
#plt.tight_layout()
#plt.savefig('/Users/pmarin/Aero_SizeDist_DOE_schematic1.png')
#
#print(gintN[ex_id,0],gcpc[ex_id,0])
#
#
#fig,ax = plt.subplots(1,1,figsize=(7,4));
#ax.plot(standDP,gdist[ex_id,:],'-k',lw=5)
#ax.plot(DP_full,gdist_extrap_cpc_save[ex_id,:],'--y',lw=3)
#ax.set_xscale('log')
#ax.set_ylabel('dN (dln $D_p$)$^{-1}$')
#ax.set_xlim([0.007,14])
#ax.set_ylim([0,3000])
#ax.set_xlabel('Particle Diameter ($D_p$, $\mu$m)')
#plt.tight_layout()
#plt.savefig('/Users/pmarin/Aero_SizeDist_DOE_schematic2.png')


fig,ax = plt.subplots(3,1,figsize=(7,10));
plt.rcParams.update({'font.size': 14})
for j in np.arange(0,3,1):   

    ex_id = np.abs(np.array(gtimes[0:cnt1])-example_times[j]).argmin()
   
    ax[j].plot(standDP,gdist[ex_id,:],c='k',lw=12,label='Original')
    ax[j].plot(DP_full,gdist_extrap_save[ex_id,:],c='0.40',lw=7,label='Step 1: Extrapolation')
    ax[j].plot(DP_full,gdist_extrap_cpc_save[ex_id,:],c='0.55',lw=6,label='Step 2: CPC Efficiency Scaling')
    ax[j].plot(DP_full,gdist_extrap_cpc_lt_save[ex_id,:],c='0.70',lw=5,label='Step 3: Systematic Bias Correction ('+str(int(med_pct_dif_both[ex_id]*100))+'%)')
    ax[j].plot(DP_full,gdist_extrap_cpc_last_save[ex_id,:],c='0.80',lw=4,label='Step 4: Adjust for Remaining Difference: '+str(int(cur_cpcint_dif[ex_id]))+' | '+str(int(cur_cpcint_pdif[ex_id]*100))+'%')
    ax[j].plot(DP_full,gdist_extrap_cpc_last_est_save[ex_id,:],c='r',lw=3,label='Step 5: CPC Efficiency Scaling')
    #ax[j].plot(DP_full,gdist_extrap_cpc_last_est_save_smooth[ex_id,:],c='m',lw=1,label='Step 5.5: Smooth')
    ax[j].set_xscale('log')
    ax[j].set_ylabel('dN (dln $D_p$)$^{-1}$')
    ax[j].set_title(strs[j]+ str(gtimes[ex_id])+savename)
    ax[j].set_title(strs[j]+ str(gtimes[ex_id]))

    if j == 2:
        ax[j].set_xlabel('Diameter ($\mu m$)')
    ax[j].legend(fontsize = 8)

plt.tight_layout()
plt.savefig(pathname+'xPUB_Final_3MCorrecttoMax_'+savename+'_Processing_Example_Monthly.png')
plt.savefig(pathname+'xPUB_Final_3MCorrecttoMax_'+savename+'_Processing_Example_Monthly.eps')
plt.close(fig)


def nearest_ind(items,pivot):
    time_diff = np.abs([date - pivot for date in items])
    return time_diff.argmin(0)

id_time = nearest_ind(gtimes,datetime(2013,5,16,15))

################################################
#####
#####  Write out processed data to netcdf files
#####
################################################
##

# QC Flag 1 --> denotes when the med_pct_dif > 50%
qcf1 = np.zeros(cnt1)
for i in np.arange(0,cnt1):
    if np.abs(med_pct_dif_both[i]) > 0.5: 
        qcf1[i] = 1

qcf2 = np.zeros(cnt1)
for i in np.arange(0,cnt1):
    if np.all(np.isnan(gdist_extrap_cpc_last_est_save[i,150:])):
        qcf2[i] = 1

qcf3 = np.zeros(cnt1)
for i in np.arange(0,cnt1):
    if gcpc[i,0] < 100 or np.isnan(gcpc[i,0]):
        qcf3[i] = 1

start_date = datetime(2009,1,1,0,0,0)
end_date = datetime(2013,12,31,0,0,0)
date_basetime = datetime(1970,1,1,0,0,0)
from netCDF4 import Dataset
import netCDF4
num_times_per_file = np.empty(0)
for tttt in np.arange(0,(end_date-start_date).days):

    date_current = start_date + timedelta(days=np.float(tttt))
    date_next = start_date + timedelta(days=tttt+1.0)
    date_current_str = date_current.strftime('%Y%m%d')
    
    s_idx = np.empty(0,int)
    s_idx_hrs = np.empty(0)
    s_idx_sec = np.empty(0)
    s_idx_base = np.empty(0)

    for ssss in np.arange(0,len(gtimes)):
        if (gtimes[ssss] >= date_current) & (gtimes[ssss] < date_next):
            s_idx = np.append(s_idx,int(ssss))
            s_idx_hrs = np.append(s_idx_hrs,(gtimes[ssss]-date_current).seconds/3600)
            s_idx_sec = np.append(s_idx_sec,(gtimes[ssss]-date_current).seconds)
            s_idx_base = np.append(s_idx_base,(date_current-date_basetime).days)   
                       
    if(len(s_idx) > 0):
        num_times_per_file = np.append(num_times_per_file,len(s_idx))
   
        dataset = Dataset(pathname+'SGP_Merged_Size_Distributions_CSA_'+date_current_str+'.nc','w',format='NETCDF4_CLASSIC')

        dim_time = dataset.createDimension('time', None)
        dim_bins_CSA = dataset.createDimension('bins_CSA',len(DP_full)-10) # only up to ~14 microns
        dim_bins_SA = dataset.createDimension('bins_SA',len(DP)-10) # only up to ~14 microns
        dim_basetime = dataset.createDimension('basetime_dim',1)

        basetime = dataset.createVariable('basetime',np.int32,('basetime_dim',))
        basetime_num = dataset.createVariable('basetime_num',np.int32,('basetime_dim',))
        time_offset = dataset.createVariable('time_offset',np.float64,('time',))

        bins_cen_CSA = dataset.createVariable('bins_cen_CSA',np.float64,('bins_CSA',))
        dN_dlnDp_CSA = dataset.createVariable('dN_dlnDp_CSA',np.float64,('time','bins_CSA'))       
        intN_CSA = dataset.createVariable('intN_CSA',np.float64,('time'))

        cpc = dataset.createVariable('cpc',np.float64,('time'))
        cpc_med_pct_dif = dataset.createVariable('cpc_med_pct_dif',np.float64,('time'))

        bins_cen_SA = dataset.createVariable('bins_cen_SA',np.float64,('bins_SA',))
        dN_dlnDp_SA = dataset.createVariable('dN_dlnDp_SA',np.float64,('time','bins_SA'))       
        intN_SA = dataset.createVariable('intN_SA',np.float64,('time'))

        qc_flag1 = dataset.createVariable('qc_flag1',np.float64,('time'))
        qc_flag2 = dataset.createVariable('qc_flag2',np.float64,('time'))
        qc_flag3 = dataset.createVariable('qc_flag3',np.float64,('time'))

        str_for_nc = netCDF4.stringtochar(np.array([date_current_str],'S8'))
        str_out = np.array([date_current_str],dtype='object')

        # Time Variables
#        basetime[:] = date_current_str
        basetime[:] = int(s_idx_base[0])
        basetime_num[:] = int(date_current_str)
        time_offset[:] = s_idx_sec

        # Adjusted Variables
        bins_cen_CSA[:] = DP_full[0:len(DP_full)-10]
        dN_dlnDp_CSA[:] = gdist_extrap_cpc_last_est_save[s_idx,0:len(DP_full)-10]
        intN_CSA[:] = gdist2_intN_est[s_idx]

        bins_cen_SA[:] = DP[0:len(DP)-10]
        dN_dlnDp_SA[:] = gdist[s_idx,0:len(DP)-10]
        intN_SA[:] = gintN[s_idx,0]

        cpc_med_pct_dif[:] = med_pct_dif_both[s_idx]
        cpc[:] = gcpc[s_idx,0]
                
        # Units
        qc_flag1[:] = qcf1[s_idx]
        qc_flag2[:] = qcf2[s_idx]
        qc_flag3[:] = qcf3[s_idx]

        date_str2 = date_current.strftime('%b %d %Y')

        basetime.units = 'days since Jan. 1 1970 000000 UTC'
        basetime_num.units = 'YYYYMMDD'
        intN_SA.units = '# per cubic cm'
        cpc_med_pct_dif.units = '%'
        cpc.units = '# per cubic cm'
        qc_flag1.units = 'unitless'
        qc_flag2.units = 'unitless'
        qc_flag3.units = 'unitless'

        dataset.close()

##########################################

# Eliminate all NaN's in the dataset
gdist2 = copy.deepcopy(gdist_extrap_cpc_est)

gdist2 = gdist2[:,0:223]
DP_full2 = DP_full[0:223]
dlogDP_full2 = dlogDP_full[0:223] # eliminate data about 14 mirons
gtimes2 = []
for i in np.arange(0,cnt1):
    if ~np.isnan(gdist2[i,10]):
        gtimes2.append(gtimes[i])
gcpc = gcpc[0:cnt1]; gcpc2 = gcpc[~np.isnan(gdist2[:,10])]
gintN = gintN[0:cnt1]; gintN2 = gintN[~np.isnan(gdist2[:,10])]
gdist2_intN2 = gdist2_intN[~np.isnan(gdist2[:,10])]
gdist2 = gdist2[~np.isnan(gdist2[:,10]),:]
#gdist2[np.isnan(gdist2)] = 0
cnt2 = len(gcpc2)
curDP = DP_full2

####################################################################################
############ Make New Variables for specific integrations of certains sizes
####################################################################################
DV_conv = np.pi/6 * np.power(DP_full2,3) # Conversion to Volume Size Distribution
SA_conv = np.pi * np.power(DP_full2,2) # Conversion to Volume Size Distribution

bin_edges = [0, 0.030, 0.140, 0.800, 30] # bin edges in microns

gintN_2_sum = np.zeros(cnt2); gintN_2 = np.zeros((len(bin_edges)-1,cnt2))
gintS_2_sum = np.zeros(cnt2); gintS_2 = np.zeros((len(bin_edges)-1,cnt2))
gintV_2_sum = np.zeros(cnt2); gintV_2 = np.zeros((len(bin_edges)-1,cnt2))
for i in np.arange(0,cnt2):
    for j in np.arange(1,len(bin_edges)):
        #curDP = copy.deepcopy(gDP[i,:])
        idx1 = np.abs(curDP-bin_edges[j-1]).argmin()
        idx2 = np.abs(curDP-bin_edges[j]).argmin()-1

        cur_dlogDP = dlogDP_full2[idx1:idx2]
        gintN_2[j-1,i] = np.nansum(gdist2[i,idx1:idx2]*cur_dlogDP) #         
        gintS_2[j-1,i] = np.nansum(gdist2[i,idx1:idx2]*SA_conv[idx1:idx2]*cur_dlogDP) #         
        gintV_2[j-1,i] = np.nansum(gdist2[i,idx1:idx2]*DV_conv[idx1:idx2]*cur_dlogDP) #         

        if j > 1 and j < 4 and (gintN_2[j-1,i] == 0 or np.isnan(gintN_2[j-1,i])): # Not for smallest size bin (since there actually could be zeros particles in this bin / especially with corrections)
            gintN_2[:,i] = np.nan;
            gintS_2[:,i] = np.nan;
            gintV_2[:,i] = np.nan;

    gintN_2_sum[i] = np.nansum(gintN_2[:,i])
    gintS_2_sum[i] = np.nansum(gintS_2[:,i])
    gintV_2_sum[i] = np.nansum(gintV_2[:,i])

gintN_2[gintN_2 < 0.0] = np.nan;
gintS_2[gintN_2 < 0.0] = np.nan;
gintV_2[gintN_2 < 0.0] = np.nan;

gintS_2[3,gintN_2[3,:] < 0.00000000001] = np.nan;
gintV_2[3,gintN_2[3,:] < 0.00000000001] = np.nan;
gintN_2[3,gintN_2[3,:] < 0.00000000001] = np.nan;

gintV_2_sum[gintN_2_sum < 10] = np.nan
gintS_2_sum[gintN_2_sum < 10] = np.nan
gintN_2_sum[gintN_2_sum < 10] = np.nan


# Number of samples
#np.count_nonzero(~np.isnan(df_120M.intN_v2))

#fig = plt.figure(figsize=(10,6));
#plt.scatter(gtimes2,gintN_2_sum,label='fin')
#plt.scatter(gtimes,gdist2_intN_est)
#plt.legend()

# For Comparison with NPF Campaiagn
bin_edges2 = [0, 0.030, 0.140] # bin edges in microns
gintN_sum_npfC = np.zeros(cnt2); gintN_npfC = np.zeros((len(bin_edges)-1,cnt2))
gintN_sum_npfU = np.zeros(cnt2); gintN_npfU = np.zeros((len(bin_edges)-1,cnt2))
gintN_npfU_E = np.zeros((len(bin_edges)-1,cnt2))
for i in np.arange(0,cnt2):
    for j in np.arange(1,len(bin_edges2)):
        #curDP = copy.deepcopy(gDP[i,:])
        idx1 = np.abs(curDP-bin_edges2[j-1]).argmin()
        idx2 = np.abs(curDP-bin_edges2[j]).argmin()
        idx3 = np.abs(DP-bin_edges2[j-1]).argmin()
        idx4 = np.abs(DP-bin_edges2[j]).argmin()

        cur_dlogDP = dlogDP_full[idx1:idx2]
        #cur_dlogDP = np.append(cur_dlogDP,cur_dlogDP[len(cur_dlogDP)-1])
        gintN_npfC[j-1,i] = np.nansum(gdist2[i,idx1:idx2]*cur_dlogDP) # 
        
        id_gtimes = np.abs(np.array(gtimes)-gtimes2[i]).argmin()


        gintN_npfU_E[j-1,i] = np.nansum(gdist_extrap_save[id_gtimes,idx1:idx2]*cur_dlogDP) #          
        gintN_npfU[j-1,i] = np.nansum(gdist[id_gtimes,idx3:idx4]*dlogDP[idx3:idx4]) #          

    gintN_sum_npfC[i] = np.nansum(gintN_npfC[:,i])
    gintN_sum_npfU[i] = np.nansum(gintN_npfU[:,i])

fig = plt.figure(figsize=(10,10))
plt.plot(gtimes2,gintN_npfU[1,:],label='Uncorrected')  
plt.plot(gtimes2,gintN_npfC[1,:],label='Corrected')   
plt.plot(gtimes,gintN,label='IntN_Orig')   
plt.plot(gtimes,gintN_full,label='IntN_Full')   

plt.legend()
###########################################################################################################################
##### Create a pandas data frame that can be used by the Power_Sectrum Code for determining spectrum for the TDMA data
###########################################################################################################################3
df_act = pd.DataFrame.from_dict({'time': gtimes2,
                                 'intN_v1':gdist2_intN2,
                                 'intN': gintN2[0:cnt2,0],
                                 'intN0': gintN_2[0,:],
                                 'intN1': gintN_2[1,:],
                                 'intN2': gintN_2[2,:],
                                 'intN3': gintN_2[3,:],
                                 'intS0': gintS_2[0,:],
                                 'intS1': gintS_2[1,:],
                                 'intS2': gintS_2[2,:],
                                 'intS3': gintS_2[3,:],
                                 'intV0': gintV_2[0,:],
                                 'intV1': gintV_2[1,:],
                                 'intV2': gintV_2[2,:],
                                 'intV3': gintV_2[3,:],
                                 'intN_v2': gintN_2_sum[:],
                                 'intS_v2': gintS_2_sum[:],
                                 'intV_v2': gintV_2_sum[:],
                                 'cpc_dma': gcpc2[0:cnt2,0],
                                 'intN0_U': gintN_npfU[0,0:cnt2],
                                 'intN1_U': gintN_npfU[1,0:cnt2],
                                 'intN0_UE': gintN_npfU_E[0,0:cnt2],
                                 'intN1_UE': gintN_npfU_E[1,0:cnt2]})
df_act = df_act.set_index('time')

hhh = df_act.index.hour; mmm = df_act.index.month; yyy = df_act.index.year;
df_act['yr'] = yyy;    df_act['mo'] = mmm;    df_act['hr'] = hhh

for i in np.arange(0,len(curDP),1):
    exec('df_act["d'+str(i)+'"] = pd.Series(gdist2[0:cnt2,'+str(i)+'], index=df_act.index)')

df_120M = df_act.resample('120Min',loffset='60Min').mean()
new_date_index_2H = [datetime(2009,1,1,1,0,0) + timedelta(hours=2*x) for x in range(0,26280)]
df_120M = df_120M.reindex(new_date_index_2H, fill_value = np.nan)
df_120M['yr'] = df_120M.index.year; df_120M['mo'] = df_120M.index.month; df_120M['hr'] = df_120M.index.hour

##### Only include this line when calculating the median distributions!
df_120M = df_120M[(df_120M.index < datetime(2011,4,14,7)) | (df_120M.index > datetime(2011,9,30,3)) ] # Eliminate period around summer 2011 with no coarse mode data

df_1440M = df_act.resample('1440Min',loffset='720Min').mean()
new_date_index_24H = [datetime(2009,1,1,12,0,0) + timedelta(hours=24*x) for x in range(0,2190)]
df_1440M = df_1440M.reindex(new_date_index_24H, fill_value = np.nan)
df_1440M['yr'] = df_1440M.index.year; df_1440M['mo'] = df_1440M.index.month; df_1440M['hr'] = df_1440M.index.hour

# Save final Pandas DataFrames for 
#pathname = '/Volumes/Data/Act_Spectrum/H/'
#pathname = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Code/H2/'


import pickle
pickle.dump([gtimes,gdist2,gintN_2,gintV_2,gintN_2_sum,gintV_2_sum,gcpc],open(pathname+"Final_Integrated_Data.p","wb"))
pickle.dump([df_120M,DP_full2],open(pathname+'df_120M_Fn'+savename+'.p',"wb"))
pickle.dump([df_1440M,DP_full2],open(pathname+'df_1440M_Fn'+savename+'.p',"wb"))

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

df_120M_DJF = pd.concat([df_120M_D,df_120M_JF])

###### Save the Median Distributions
###### Be Careful.. Realize whether need to drop the summer 2011 data
#dn = np.zeros((len(DP_full2),5))
#ds = np.zeros((len(DP_full2),5))
#dv = np.zeros((len(DP_full2),5))
#
#dn[:,0] = np.nanpercentile(df_120M.iloc[:,bin0_id:],50,axis=0)
#dn[:,1] = np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:],50,axis=0)
#dn[:,2] = np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:],50,axis=0)
#dn[:,3] = np.nanpercentile(df_120M_SON.iloc[:,bin0_id:],50,axis=0)
#dn[:,4] = np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:],50,axis=0)
#
#ds[:,0] = np.nanpercentile(df_120M.iloc[:,bin0_id:]*SA_conv,50,axis=0)
#ds[:,1] = np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*SA_conv,50,axis=0)
#ds[:,2] = np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*SA_conv,50,axis=0)
#ds[:,3] = np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*SA_conv,50,axis=0)
#ds[:,4] = np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*SA_conv,50,axis=0)
#  
#dv[:,0] = np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,50,axis=0)
#dv[:,1] = np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*DV_conv,50,axis=0)
#dv[:,2] = np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*DV_conv,50,axis=0)
#dv[:,3] = np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*DV_conv,50,axis=0)
#dv[:,4] = np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*DV_conv,50,axis=0)
#
#pickle.dump([DP_full2,dn,ds,dv],open(pathname+"Median_Size_Distributions.p","wb"))

#######################################################################################
#######################################################################################
###########
###########                        PLOT TIME
###########
#######################################################################################
#######################################################################################

#mins = np.zeros(len(AAA))
#AAA = np.diff(gtimes)
#for i in np.arange(0,len(mins)):
#    mins[i] = AAA[i].seconds/60
#    
#    
#bins=np.arange(0,100,2)
#np.nanmedian(mins)
#fig = plt.figure(figsize=(10,5))
#plt.hist(mins,bins=bins)
#np.nanpercentile(mins,25)


#savepath = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Code/H2'
savepath = pathname

######## Plot all aerosol data on a yearly basis ( Paper Quality)
#import matplotlib   
#import matplotlib.dates as mdates
#
#plt.rcParams.update({'font.size': 12})
#fig,ax = plt.subplots(5,1,figsize=(12,13))
#cnt_plt = 0
#clabels = np.arange(-5.0,5.1,0.1)
#bin0_id = 25
#for t in np.arange(2009,2014,1):
#    
#    xlims = matplotlib.dates.date2num([datetime(t,1,1), datetime(t,12,31)])
#
#    #cf = ax[cnt_plt].contourf(gtimes,DP,np.log10(np.transpose(gdist2[0:len(gtimes),:])),50,lw='none',cmap=plt.cm.jet)
#    cf = ax[cnt_plt].contourf(df_120M.index,DP_full2,np.log10(np.transpose(df_120M.iloc[:,bin0_id:])),levels=clabels,lw='none', extend='both', cmap=plt.cm.nipy_spectral)
#    ax[cnt_plt].set_ylabel('Diameter ($\mu$m)')
#    ax[cnt_plt].set_yscale('log')
#    ax[cnt_plt].set_xlim(xlims)
#    ax[cnt_plt].set_ylim([0.007,13])
#    cb = plt.colorbar(cf,ax=ax[cnt_plt],pad=0.1,ticks=[-4,-3,-2,-1,0,1,2,3,4])
#    cb.ax.set_yticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','1','10','$10^{2}$','$10^{3}$','$10^{4}$'])
#    cb.ax.set_ylabel('dN dln$(D_{p})^{-1}$')
#    ax[cnt_plt].set_title(str(t),loc='left')
#
#    #ax.set_xticks(xdates)
#    ax[cnt_plt].xaxis.set_major_formatter(mdates.DateFormatter("%m|%y"))
#
#    ax2 = ax[cnt_plt].twinx()
#    ax2.plot(df_120M.index,df_120M.intN_v2,'ok',ms='0.5',label='DMA|CPC_v2')
#    ax2.set_xlim(xlims)
#    ax2.set_ylim([0,20000])
#    ax2.set_ylabel('Number Conc. (# cm$^{-3}$)')
#
#
#    cnt_plt = cnt_plt + 1;
#plt.tight_layout()
#plt.savefig(savepath+'FINAL_Corrected_Aerosol_Temporal_Data_Summary_Narrow'+savename+'.png')
#plt.savefig(savepath+'FINAL_Corrected_Aerosol_Temporal_Data_Summary_Narrow'+savename+'.eps')
##plt.savefig(savepath+'FINAL_Corrected_Aerosol_Temporal_Data_Summary_Narrow'+savename+'.pdf')
#

####### Plot mean aerosol size distribution by season ( Paper Quality)
import matplotlib   
import matplotlib.dates as mdates

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color=['black','blue','crimson','goldenrod','0.75'])

bin0_id = 25
pcttile = 50

plt.rcParams.update({'font.size': 17})
fig,ax = plt.subplots(4,1,figsize=(7,14))
# Specify bin edges
max_val = 3500
#ax[0].plot([bin_edges[1],bin_edges[1]],[0,max_val],c='0.6',lw=5)
#ax[0].plot([bin_edges[2],bin_edges[2]],[0,max_val],c='0.6',lw=5)
#ax[0].plot([bin_edges[3],bin_edges[3]],[0,max_val],c='0.6',lw=5)
ax[0].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:],pcttile,axis=0),'-',lw=7,label='ALL')
ax[0].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:],pcttile,axis=0),'-',lw=4,label='MAM')
ax[0].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:],pcttile,axis=0),'-',lw=4,label='JJA')
ax[0].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:],pcttile,axis=0),'-',lw=4,label='SON')
ax[0].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:],pcttile,axis=0),'-',lw=4,label='DJF')

np.nansum(np.nanpercentile(df_120M.iloc[:,bin0_id:],pcttile,axis=0)*dlogDP_full2)
np.nanmedian(gcpc2)
np.nanmedian(gdist2_intN_est)
np.nanmedian(gintN_2_sum)


pcttile_l = 25
pcttile_h = 75
ax[0].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:],pcttile_l,axis=0),':',lw=3.5)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:],pcttile_l,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:],pcttile_l,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:],pcttile_l,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:],pcttile_l,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:],pcttile_h,axis=0),':',lw=3.5)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:],pcttile_h,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:],pcttile_h,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:],pcttile_h,axis=0),':',lw=2)
ax[0].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:],pcttile_h,axis=0),':',lw=2)

#ax[0,0].set_xscale('log')
#ax[0,0].set_xlabel('Diameter ($\mu$m)')
ax[0].set_ylabel('dN dlnD$_{p}^{-1}$ \n (# cm$^{-3}$)')
ax[0].grid()
#ax[0,0].set_yscale('log')
ax[0].set_xscale('log')
#ax[0].set_yscale('log')
ax[0].set_ylim([0.001,max_val])
ax[0].set_xlim([.007,12.5])
ax[0].legend(fontsize=12,loc="upper right")
ax[0].text(0.015,max_val*0.85,'(a)')

max_val2 = 120
#ax[1].plot([bin_edges[1],bin_edges[1]],[0,max_val2],c='0.6',lw=5)
#ax[1].plot([bin_edges[2],bin_edges[2]],[0,max_val2],c='0.6',lw=5)
#ax[1].plot([bin_edges[3],bin_edges[3]],[0,max_val2],c='0.6',lw=5)
ax[1].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:]*SA_conv,pcttile,axis=0),'-',lw=7,label='ALL')
ax[1].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*SA_conv,pcttile,axis=0),'-',lw=4,label='MAM')
ax[1].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*SA_conv,pcttile,axis=0),'-',lw=4,label='JJA')
ax[1].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*SA_conv,pcttile,axis=0),'-',lw=4,label='SON')
ax[1].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*SA_conv,pcttile,axis=0),'-',lw=4,label='DJF')

pcttile_l = 25
pcttile_h = 75
ax[1].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:]*SA_conv,pcttile_l,axis=0),':',lw=3.5)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*SA_conv,pcttile_l,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*SA_conv,pcttile_l,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*SA_conv,pcttile_l,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*SA_conv,pcttile_l,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:]*SA_conv,pcttile_h,axis=0),':',lw=3.5)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*SA_conv,pcttile_h,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*SA_conv,pcttile_h,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*SA_conv,pcttile_h,axis=0),':',lw=2)
ax[1].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*SA_conv,pcttile_h,axis=0),':',lw=2)

ax[1].set_ylabel('dS dlnD$_{p}^{-1}$ \n ($\mu$m$^{2}$ cm$^{-3}$)')
ax[1].set_ylim([0.001,max_val2])
ax[1].set_xlim([.007,12.5])
ax[1].grid()
ax[1].set_xscale('log')
#ax[1].set_yscale('log')
ax[1].text(0.015,max_val2*0.85,'(b)')

max_val2 = 5
#ax[2].plot([bin_edges[1],bin_edges[1]],[0,max_val2],c='0.6',lw=5)
#ax[2].plot([bin_edges[2],bin_edges[2]],[0,max_val2],c='0.6',lw=5)
#ax[2].plot([bin_edges[3],bin_edges[3]],[0,max_val2],c='0.6',lw=5)
DV_conv = np.pi/6 * np.power(DP_full2,3)
ax[2].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0),'-',lw=7,label='ALL')
ax[2].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0),'-',lw=4,label='MAM')
ax[2].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0),'-',lw=4,label='JJA')
ax[2].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0),'-',lw=4,label='SON')
ax[2].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0),'-',lw=4,label='DJF')

pcttile_l = 25
pcttile_h = 75
ax[2].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile_l,axis=0),':',lw=3.5)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*DV_conv,pcttile_l,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*DV_conv,pcttile_l,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*DV_conv,pcttile_l,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*DV_conv,pcttile_l,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile_h,axis=0),':',lw=3.5)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*DV_conv,pcttile_h,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*DV_conv,pcttile_h,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*DV_conv,pcttile_h,axis=0),':',lw=2)
ax[2].plot(DP_full2,np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*DV_conv,pcttile_h,axis=0),':',lw=2)

ax[2].set_ylabel('dV dlnD$_{p}^{-1}$ \n ($\mu$m$^{3}$ cm$^{-3}$)')
ax[2].set_ylim([0.001,max_val2])
ax[2].set_xlim([.007,12.5])
ax[2].grid()
ax[2].set_xscale('log')
#ax[1].set_yscale('log')
ax[2].text(0.015,max_val2*0.85,'(c)')

max_val2 =8
#ax[3].plot([bin_edges[1],bin_edges[1]],[-60,60],c='0.6',lw=5)
#ax[3].plot([bin_edges[2],bin_edges[2]],[-60,60],c='0.6',lw=5)
#ax[3].plot([bin_edges[3],bin_edges[3]],[-60,60],c='0.6',lw=5)
ax[3].plot(DP_full2,(np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)-np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0))/np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)*100,'-',lw=6,label='All')
ax[3].plot(DP_full2,(np.nanpercentile(df_120M_MAM.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)-np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0))/np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)*100,'-',lw=4,label='MAM')
ax[3].plot(DP_full2,(np.nanpercentile(df_120M_JJA.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)-np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0))/np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)*100,'-',lw=4,label='JJA')
ax[3].plot(DP_full2,(np.nanpercentile(df_120M_SON.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)-np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0))/np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)*100,'-',lw=4,label='SON')
ax[3].plot(DP_full2,(np.nanpercentile(df_120M_DJF.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)-np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0))/np.nanpercentile(df_120M.iloc[:,bin0_id:]*DV_conv,pcttile,axis=0)*100,'-',lw=4,label='DJF')
ax[3].set_xscale('log')
ax[3].set_ylabel('Percentage Difference \n in Median Values')
ax[3].set_ylim([-60,60])
ax[3].set_xlim([.007,12.5])
ax[3].grid()
ax[3].set_xlabel('Diameter ($\mu$m)')
ax[3].text(0.015,43,'(d)')

plt.tight_layout()

plt.savefig(savepath+'Final_e4'+savename+'_2W_vCPC_XPUB_Corrected_Log_'+str(pcttile)+'th_Aerosol_Size_Distributions_Taller_wDV_nosummer.png')
plt.savefig(savepath+'Final_e4'+savename+'_2W_vCPC_XPUB_Corrected_Log_'+str(pcttile)+'th_Aerosol_Size_Distributions_Taller_wDV_nosummer.eps')

#plt.savefig('/Users/pmarin/Aerosol_Size_Distributions_Taller_wDV_NoLines.png')


#plt.savefig(savepath+'Final_e4'+savename+'_2W_vCPC_XPUB_Corrected_Log_'+str(pcttile)+'th_Aerosol_Size_Distributions_Taller_wDV_Drop.png')
#plt.savefig(savepath+'Final_e4'+savename+'_2W_vCPC_XPUB_Corrected_Log_'+str(pcttile)+'th_Aerosol_Size_Distributions_Taller_wDV_Drop.eps')

####################################################################################
############### Compare data with NPF Field Campaign Data
############### Need to also run the Read_NPF code to get the dataframe dfn_120M
############### Run /Volumes/Data/Act_Spectrum/FFFF/Read_NPF_data_vFinal.py
###################################################################################
import matplotlib.dates as mdates
plt.rcParams.update({'font.size': 11})

fig,ax = plt.subplots(2,1,figsize=(15,7))
ylims = [50,50000]

xlims = [datetime(2013,4,19,12,0,0), datetime(2013,5,3,12,0,0)]
ax[0].plot(dfn_120M.index,dfn_120M['nm730'],'-',c='crimson',lw=6,label='NPFS')
ax[0].plot(df_120M.index,df_120M.intN0_U,'-k',lw=3,label='TDMA-Original')
#ax[0].plot(df_120M.index,df_120M.intN0_UE,':k',lw=5,label='TDMA-Extrap')
ax[0].plot(df_120M.index,df_120M.intN0,'-',c='0.50',lw=3,label='TDMA-Adjusted')
ax[0].set_ylim(ylims);
ax[0].set_ylim(ylims); ax[0].set_yscale('log')
ax[0].set_ylabel('Number Concentration \n (# $cm^{-3}$)')

xdates = [xlims[0] + timedelta(hours=x) for x in range(0,1000,12)]
ax[0].set_xticks(xdates)
ax[0].set_xlim(xlims)
ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
ax[0].set_xlabel('UTC Time (Day|Hour)')
ax[0].legend(loc="bottom left",fontsize=10)
ax[0].grid()

xlims = [datetime(2013,5,3,12,0,0), datetime(2013,5,17,12,0,0)]
ax[1].plot(dfn_120M.index,dfn_120M['nm730'],'-',c='crimson',lw=6,label='NPFS')
ax[1].plot(df_120M.index,df_120M.intN0_U,'-k',lw=3,label='TDMA-Original')
#ax[1].plot(df_120M.index,df_120M.intN0_UE,':k',lw=5,label='TDMA-Extrap')
ax[1].plot(df_120M.index,df_120M.intN0,'-',c='0.50',lw=3,label='TDMA-Adjusted')
ax[1].set_ylim(ylims);
ax[1].set_ylim(ylims); ax[1].set_yscale('log')
xdates = [xlims[0] + timedelta(hours=x) for x in range(0,1000,12)]
ax[1].set_xticks(xdates)
ax[1].set_ylabel('Number Concentration \n (# $cm^{-3}$)')
ax[1].set_xlim(xlims)
ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
#ax[1].legend(loc="bottom left")
ax[1].set_xlabel('UTC Time (Day|Hour)')
ax[1].grid()
plt.tight_layout()

plt.savefig(pathname+'Comparison_with_NPF2013Data_730_3MCorrecttoMax_'+savename+'_curvefitDP2_meandiff.png')
plt.savefig(pathname+'Comparison_with_NPF2013Data_730_3MCorrecttoMax_'+savename+'_curvefitDP2_meandiff.eps')
#plt.close(fig)

fig,ax = plt.subplots(2,1,figsize=(15,7))
ylims = [100,10000]

xlims = [datetime(2013,4,19,12,0,0), datetime(2013,5,3,12,0,0)]
ax[0].plot(dfn_120M.index,dfn_120M['nm30140'],'-',c='crimson',lw=9,label='NPFS')
ax[0].plot(df_120M.index,df_120M.intN1_U,'-k',lw=5,label='TDMA-Original')
ax[0].plot(df_120M.index,df_120M.intN1,c='0.50',lw=5,label='TDMA-Adjusted')
ax[0].set_ylim(ylims);
ax[0].set_ylim(ylims); ax[0].set_yscale('log')

xdates = [xlims[0] + timedelta(hours=x) for x in range(0,1000,12)]
ax[0].set_xticks(xdates)
ax[0].set_ylabel('(# $cm^{-3}$)')
ax[0].set_xlim(xlims)
ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
ax[0].legend(loc="bottom left",fontsize=14)
ax[0].set_xlabel('UTC Time (Day|Hour)')
ax[0].grid()

xlims = [datetime(2013,5,3,12,0,0), datetime(2013,5,17,12,0,0)]
ax[1].plot(dfn_120M.index,dfn_120M['nm30140'],'-',c='crimson',lw=8,label='NPFS')
ax[1].plot(df_120M.index,df_120M.intN1_U,'-k',lw=5,label='TDMA-Original')
ax[1].plot(df_120M.index,df_120M.intN1,c='0.50',lw=5,label='TDMA-Adjusted')
ax[1].set_ylim(ylims);
ax[1].set_ylim(ylims); ax[1].set_yscale('log')

xdates = [xlims[0] + timedelta(hours=x) for x in range(0,1000,12)]
ax[1].set_xticks(xdates)
ax[1].set_ylabel('(# $cm^{-3}$)')
ax[1].set_xlim(xlims)
ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%d|%H"))
#ax[1].legend(loc="bottom left")
ax[1].set_xlabel('UTC Time (Day|Hour)')
ax[1].grid()
plt.tight_layout()

plt.savefig(pathname+'Comparison_with_NPF2013Data_30140_3MCorrecttoMax_'+savename+'_curvefitDP2_meandiff.png')
plt.savefig(pathname+'Comparison_with_NPF2013Data_30140_3MCorrecttoMax_'+savename+'_curvefitDP2_meandiff.eps')
#plt.close(fig)

dt1 = datetime(2013,4,19,1,0,0)
dt2 = datetime(2013,5,23,23,0,0)

id_npf1 = np.abs(dfn_120M.index-dt1).argmin()
id_npf2 = np.abs(dfn_120M.index-dt2).argmin()

id_dma1 = np.abs(df_120M.index-dt1).argmin()
id_dma2 = np.abs(df_120M.index-dt2).argmin()

aa1 = df_120M.intN0_U[id_dma1:id_dma2]
aa2 = df_120M.intN0[id_dma1:id_dma2]
aa3 = df_120M.intN1_U[id_dma1:id_dma2]
aa4 = df_120M.intN1[id_dma1:id_dma2]

bb1 = dfn_120M.nm730[id_npf1:id_npf2]
bb2 = dfn_120M.nm30140[id_npf1:id_npf2]

###########################################################################################################################3
df_corr = pd.DataFrame.from_dict({'time': aa1.index,
                                 'intN0_U':aa1.values,
                                 'intN1_U':aa3.values,
                                 'intN0':aa2.values,
                                 'intN1':aa4.values,
                                 'intN0_npf':bb1.values,
                                 'intN1_npf':bb2.values})
df_corr = df_corr.set_index('time')

df_corr2 = pd.DataFrame.from_dict({'time': aa1.index,
                                 'intN0_U':standardize(aa1.values),
                                 'intN1_U':standardize(aa3.values),
                                 'intN0':standardize(aa2.values),
                                 'intN1':standardize(aa4.values),
                                 'intN0_npf':standardize(bb1.values),
                                 'intN1_npf':standardize(bb2.values)})
df_corr2 = df_corr2.set_index('time')

#fig = plt.figure(figsize=(7,7));
#plt.scatter(aa3,bb2,label='npfU')
#plt.scatter(aa4,bb2,label='npfC')
#plt.legend()
#plt.xlim([0,14000])
#plt.ylim([0,14000])
#plt.xlabel('NPF Data')
#plt.ylabel('DMA Data')
#plt.tight_layout()

# Get Correlation values for these variables
df_corr.corr()
df_corr2.corr()