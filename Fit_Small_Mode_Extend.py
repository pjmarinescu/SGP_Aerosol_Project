#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:19:10 2018
@author: pmarin
"""

import matplotlib
matplotlib.use("QT4agg")

# Calculate Power Spectrum of Data
import numpy as np
import copy
import matplotlib.pyplot as plt
import datetime
from netCDF4 import Dataset
import glob
from scipy import optimize
from scipy.signal import argrelextrema

# Get dimensions of the data
plot_on = 0
save_on = 1
pathname = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/BOTH_CPC_CCN_FINAL/Data/Data_1nm/'
plotsave = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/BOTH_CPC_CCN_FINAL/Data/Plots_v2/'
savedir = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/BOTH_CPC_CCN_FINAL/Data/Nucl_Mode_LogN_Files_Final/'

files = glob.glob(pathname+'*.nc')

parN = []
parS = []
parD = []

for f in np.arange(4,len(files)):
#for f in np.arange(1274,1275):
#for f in np.arange(356,357):
    #for f in np.arange(0,1):        

    print('File Number:'+str(f))        
    # Read in CSA size distribution data
    file_now = files[f]
    fi = Dataset(file_now, "a", format="NETCDF4")
    dNdlnDP = fi['dN_dlnDp_CSA'][:]
    qc1 = fi['qc_flag1'][:]
    qc2 = fi['qc_flag2'][:]
    qc3 = fi['qc_flag3'][:]
    print(file_now)
    
    Dp = fi['bins_cen_CSA'][:]

    dfsgs

    # Create extended size range down to 1 nm
    dp_n = copy.deepcopy(Dp)
    for d in np.arange(0,100):       
        log_diff = np.log(dp_n[1])-np.log(dp_n[0])
        new_dp = np.exp(np.log(dp_n[0])-log_diff)
        if new_dp < 0.001:
            break
        else:
            dp_n = np.insert(dp_n,0,new_dp)            
    
    # Calculate bin edges for calculating dlnDP
    dpn_edges = (dp_n[0:len(dp_n)-1]+dp_n[1:])/2
    dpn_edges = np.concatenate((dpn_edges[0]-(np.diff(dpn_edges[0:2])),dpn_edges,dpn_edges[-1]+(np.diff(dpn_edges[-2:]))),axis=0)
    dlnDP = np.diff(np.log(dpn_edges))

    base_time = fi['basetime'][:]
    current_date = datetime.datetime(1970,1,1,0,0,0) + datetime.timedelta(days=int(base_time[0]))
    time_offset = fi['time_offset'][:]
    qc3 = fi['qc_flag3'][:]   
    intN = np.zeros(len(time_offset)); intN[:] = np.nan
    dNdlnDp_f = np.zeros((len(time_offset),len(dp_n))); dNdlnDp_f[:] = np.nan
    
#    sadfad
    for t in np.arange(0,len(time_offset)):  

        if qc3[t] == 1:
            print('No CPC data; No CSA size distribution')
            continue

        else:
            time_now = current_date + datetime.timedelta(seconds=int(time_offset[t]))
            print(time_now)
            
            ## create string from date for savename
            yy = str(time_now.year)
            mm = str(time_now.month)
            if(len(mm)) == 1:
                mm = '0'+mm
            dd = str(time_now.day)       
            if(len(dd)) == 1:
                dd = '0'+dd
            ho = str(time_now.hour)       
            if(len(ho)) == 1:
                ho = '0'+ho
            mi = str(time_now.minute)       
            if(len(mi)) == 1:
                mi = '0'+mi
            se = str(time_now.second)       
            if(len(se)) == 1:
                se = '0'+se

            savename = 'FSD_'+yy+mm+dd+'_'+ho+mi+se+'_LogN.png'
    
    # For testing
    #        if ho+mi+se == '210712':
    #            adfafasdf
    
            # Find part of initial aerosol size distribution to use for fitting
            # Try to use local minimum that is not at the very left edge of the size distribution        
            min_ids = argrelextrema(dNdlnDP[t,:],np.less)
            if len(min_ids[0]) == 0: # if no local minimum, use 20 nm
                id_dp = np.argmin(np.abs(Dp-0.015))
            else: # Create smaller buffer to the left of local minimum, because using the actuall local minimum often doesn't provide the greatest fit (not steep enough)
                if (min_ids[0][0] < 5) and (len(min_ids[0]) > 1):
                    id_dp = min_ids[0][1]
                elif (min_ids[0][0] > 25):
                    id_dp = min_ids[0][0] - 15
                elif (min_ids[0][0] > 20):
                    id_dp = min_ids[0][0] - 10
                elif (min_ids[0][0] > 16):
                    id_dp = min_ids[0][0] - 8
                elif (min_ids[0][0] > 12):
                    id_dp = min_ids[0][0] - 5
                elif (min_ids[0][0] > 9):
                    id_dp = min_ids[0][0] - 3
                else:
                    id_dp = min_ids[0][0]
    
            # Don't allow the use of data larger than 20 nm (getting out of the nucleation mode)
            if id_dp > 24:
                id_dp = 23
            if id_dp < 5:
                id_dp = np.argmin(np.abs(Dp-0.015))
    
            dNdlnDp_s = dNdlnDP[t,0:id_dp]
            dp_s = Dp[0:id_dp]
            dlogdp_s = np.diff(np.log(Dp[0:id_dp+1]))
            
            #DP_full_mids = (DP_full[0:len(DP_full)-1]+DP_full[1:])/2
#DP_full_mids = np.concatenate((DP_full_mids[0]-(np.diff(DP_full_mids[0:2])),DP_full_mids,DP_full_mids[-1]+(np.diff(DP_full_mids[-2:]))),axis=0)
#dlogDP_full = np.diff(np.log(DP_full_mids))

            
            
            # If there are no particles at the smallest end of the size distribution, skip fitting analysis
            intN_s = np.nansum(dNdlnDp_s[0:2]*dlogdp_s[0:2])
            if intN_s > 0:
    
                if ( dNdlnDp_s[0] == np.nanmin(dNdlnDp_s[0:]) ) and ( dNdlnDp_s[0] < (0.5 *np.nanmean(dNdlnDp_s[0:])) )  :
    
                    print(dNdlnDp_s[0],np.nanmin(dNdlnDp_s))
                    # Fit polynomial
                    dNdlnDp_s_now = dNdlnDP[t,0:11]
                    dp_s_now = Dp[0:11]
    
    #                fitfunc = lambda p, d: (p[0] + (p[1] * d) + (p[2] * np.power(d,2.0)))  # Target function
    #                fitfunc = lambda p, d: (p[0] * np.exp(p[1]*d))  # Target function
    #                errfunc = lambda p, d, y: fitfunc(p, d) - y # Distance to the target function
    #                p0 = [10,-2]
    #                out = optimize.least_squares(errfunc, p0[:], args=(dp_s_now,dNdlnDp_s_now))
    #                print(out.x)
    
                    fitfunc = lambda p, d: (p[0] / (np.log(p[1])*np.sqrt(2*np.pi))) * np.exp(-(np.power(np.log(d)-np.log(p[2]),2.0))/(2*np.power(np.log(p[1]),2.0))) # Target function
                    errfunc = lambda p, d, y: fitfunc(p, d) - y # Distance to the target function
                    p0 = [3000, 2.5, 0.006] # Initial guess for the parameters
                    bounds = ([0,1.0,0.004],[1E5,3.0,0.015])
                    out = optimize.least_squares(errfunc, p0[:], bounds=bounds, args=(dp_s_now,dNdlnDp_s_now))
                    print(out.x)
    
                    dNdlnDp_f[t,0:55] = fitfunc(out.x,dp_n)[0:55]
                    dNdlnDp_f[t,55:] = dNdlnDP[t,:]
    
                else:                            
                    # Fit the first set
                    fitfunc = lambda p, d: (p[0] / (np.log(p[1])*np.sqrt(2*np.pi))) * np.exp(-(np.power(np.log(d)-np.log(p[2]),2.0))/(2*np.power(np.log(p[1]),2.0))) # Target function
                    errfunc = lambda p, d, y: fitfunc(p, d) - y # Distance to the target function
                    p0 = [3000, 2.5, 0.006] # Initial guess for the parameters
                    bounds = ([0,1.0,0.004],[1E5,3.0,0.015])
                    out = optimize.least_squares(errfunc, p0[:], bounds=bounds, args=(dp_s,dNdlnDp_s))
                    print(out.x)
                    
                    dNdlnDp_f[t,0:55] = fitfunc(out.x,dp_n)[0:55]
                    dNdlnDp_f[t,55:] = dNdlnDP[t,:]
                
                if plot_on == 1:
                    fig = plt.figure(figsize=(7,5))
                    plt.plot(Dp,dNdlnDP[t,:],'y-',lw=6,label='CSA Size Dist.')  
                    if ( dNdlnDp_s[0] == np.nanmin(dNdlnDp_s[0:]) ) and ( dNdlnDp_s[0] < (0.5 *np.nanmean(dNdlnDp_s[0:])) )  :
                        plt.plot(dp_n,fitfunc(out.x,dp_n),'r-',lw=4,label='LogN Fit -- Nt:'+str(int(out.x[0]))+' cm$^{-3}$, sig:'+str(np.round(out.x[1],2))+', Dpm:'+str(np.round(out.x[2]*1000,1))+' nm')
        #                plt.plot(dp_n,fitfunc(out.x,dp_n),'r-',lw=4,label=str(out.x[0])+'* exp('+str(out.x[1])+')')
        #                plt.plot(dp_n,fitfunc(out.x,dp_n),'r-',lw=4,label=str(int(out.x[0]))+'+ '+str(int(out.x[1]))+'*Dp + '+str(int(out.x[2]))+'*Dp$^{2}$')
                        plt.plot(dp_s_now,dNdlnDp_s_now,'cx',label='CSA Size Dist. Used for Fit')
                    else:
                        plt.plot(dp_n,fitfunc(out.x,dp_n),'r-',lw=4,label='LogN Fit -- Nt:'+str(int(out.x[0]))+' cm$^{-3}$, sig:'+str(np.round(out.x[1],2))+', Dpm:'+str(np.round(out.x[2]*1000,1))+' nm')
                        plt.plot(dp_s,dNdlnDp_s,'cx',label='CSA Size Dist. Used for Fit')
                    plt.plot(dp_n,dNdlnDp_f,'k-',lw=2,label='Final Size Distribution')
                    plt.xscale('log')
                    plt.xlabel('Dp ($\mu$m)')
                    plt.ylabel('Dn/DlnDP')
                    plt.ylim([0,np.nanmax(dNdlnDp_f)*1.1])
                    plt.legend()
                    plt.title('idx: '+str(t))
                    #plt.savefig(plotsave+savename)
                    #plt.close(fig)
                
#                parN.append(out.x[0])
#                parS.append(out.x[1])
#                parD.append(out.x[2]*1000)
     
            else:
                print('No Particles in Small Mode')
                
                dNdlnDp_f[t,0:55] = 0
                dNdlnDp_f[t,55:] = dNdlnDP[t,:]
#    
#                fig = plt.figure(figsize=(7,5))
#                plt.plot(Dp,dNdlnDP[t,:],'y-',lw=6,label='CSA Size Dist.')  
#                #plt.plot(dp_s,dNdlnDp_s,'cx',label='CSA Size Dist. Used for Fit')
#                plt.plot(dp_n,dNdlnDp_f,'k-',lw=2,label='Final Size Distribution')
#                plt.xscale('log')
#                plt.xlabel('Dp ($\mu$m)')
#                plt.ylabel('Dn/DlnDP')
#                plt.ylim([0,np.nanmax(dNdlnDp_f)*1.1])
#                plt.legend()
#                plt.savefig(plotsave+savename)
#                plt.close(fig)
    
                #parN.append(np.nan)
                #parS.append(np.nan)
                #parD.append(np.nan)
                
            intN[t] = np.nansum(dNdlnDp_f*dlnDP)

    if save_on == 1:                    
        # Create variables        
        dim_bins_CCSA = fi.createDimension('bins_CCSA',len(dp_n)) 
        bins_cen_CCSA = fi.createVariable('bins_cen_CCSA',np.float64,('bins_CCSA',))
        intN_CCSA = fi.createVariable('intN_CCSA',np.float64,('time',))
        dN_dlnDp_CCSA = fi.createVariable('dN_dlnDp_CCSA',np.float64,('time','bins_CCSA'))   
    
        # Write Variables
        bins_cen_CCSA[:] = dp_n[:]
        dN_dlnDp_CCSA[:] = dNdlnDp_f[:,:]
        intN_CCSA[:] = intN[:]
        
        bins_cen_CCSA.units = 'microns'
        dN_dlnDp_CCSA.units = '# per cubic cm'
        intN_CCSA.units = '# per cubic cm'
        
    fi.close()
    
    
    
    
####################
### Testing Reading in Files that were adjusted in the above code
####################

##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Thu Feb  1 10:19:10 2018
#@author: pmarin
#"""
#
#import matplotlib
#matplotlib.use("QT4agg")
#
## Calculate Power Spectrum of Data
#import numpy as np
#import copy
#import matplotlib.pyplot as plt
#import datetime
#from netCDF4 import Dataset
#import glob
#from scipy import optimize
#from scipy.signal import argrelextrema
#
## Get dimensions of the data
#plot_on = 0
#pathname = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/BOTH_CPC_CCN_FINAL/Data/Data_1nm/'
#plotsave = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/BOTH_CPC_CCN_FINAL/Data/Plots_v2/'
#savedir = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/Data/BOTH_CPC_CCN_FINAL/Data/Nucl_Mode_LogN_Files_Final/'
#
#files = glob.glob(pathname+'*.nc')
#
#parN = []
#parS = []
#parD = []
#
#for f in np.arange(1,1603):
#
#    print('File Number:'+str(f))        
#    # Read in CSA size distribution data
#    file_now = files[f]
#    fi = Dataset(file_now, "a", format="NETCDF4")
#    print(file_now)
#    
#    Dp = fi['bins_cen_CCSA'][:]
#    intN = fi['intN_CCSA'][:]
#    dN = fi['dN_dlnDp_CCSA'][:]
#
#    Dp2 = fi['bins_cen_CSA'][:]
#    intN2 = fi['intN_CSA'][:]
#    dN2 = fi['dN_dlnDp_CSA'][:]
#
#    Dp3 = fi['bins_cen_SA'][:]
#    intN3 = fi['intN_SA'][:]
#    dN3 = fi['dN_dlnDp_SA'][:]
#
#    
#    time_offset = fi['time_offset'][:]
#    
#    for t in np.arange(0,len(time_offset),4):
#                
#        fig = plt.figure(figsize=(7,5))
#        plt.plot(Dp,dN[t,:],'y-',lw=6,label='CCSA Size Dist.')  
#        plt.plot(Dp2,dN2[t,:],'r-',lw=3,label='CSA Size Dist.')  
#        plt.plot(Dp3,dN3[t,:],'k-',lw=1,label='SA Size Dist.')  
#        plt.xscale('log')
#        plt.xlabel('Dp ($\mu$m)')
#        plt.ylabel('Dn/DlnDP')
#        plt.ylim([0,np.nanmax(dN[t,:])*1.1])
#        plt.title(file_now+' | '+str(t))
#        plt.legend()
#        
#    
#    fi.close()
