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
import itertools
import netCDF4 as nc
from netCDF4 import Dataset
import copy

######################################################################
# Read in PBL Height Files from ARM-Archive
######################################################################
pathname3 = '/Volumes/Data/Act_Spectrum/PBL/'  # May 2007 - Dec. 2013
#pathname3 = '/Volumes/Data/Act_Spectrum/CPC_Other/'
#g = MFDataset(pathname2+"sgpmet*cdf")
g = os.listdir(pathname3)
p_time0 = list(); pbl_ll0 = list(); pbl_r250 = list(); pres0 = list(); 
temp0 = list(); pbl_time0 = list(); qc_pbl_ll0 = list(); qc_pbl_r250 = list()
pbl_r50 = list(); qc_pbl_r50 = list()
pbl_hef0 = list(); qc_pbl_hef0 = list()
pbl_reg0 = list()
for i in g:
    print('Reading Filename:',i)
    if i.startswith('sgppblhtsonde1mcfarlC1.s1.2009') or i.startswith('sgppblhtsonde1mcfarlC1.s1.2010') or i.startswith('sgppblhtsonde1mcfarlC1.s1.2011') or i.startswith('sgppblhtsonde1mcfarlC1.s1.2012') or i.startswith('sgppblhtsonde1mcfarlC1.s1.2013'):
    #if i.startswith('sgpnoaaaosC1.b1.20090505'):
    #if i.startswith('sgpnoaaaosC1.b1.20121125'):
              
        fi = Dataset(pathname3+i, "r", format="NETCDF4")  
        #nc2 = netcdf.netcdf_file(pathname3+i,'r')
        base = nc.num2date(fi.variables['base_time'][:], units=fi.variables['base_time'].units)
        to_temp = fi.variables['time_offset'][:]    
        arr = np.array([base + timedelta(seconds=i) for i in to_temp])
        for item in arr:
            p_time0.append(item)    
        pbl_ll0.append(fi.variables['pbl_height_liu_liang'][:])
        qc_pbl_ll0.append(fi.variables['qc_pbl_height_liu_liang'][:])

        pbl_r250.append(fi.variables['pbl_height_bulk_richardson_pt25'][:])
        qc_pbl_r250.append(fi.variables['qc_pbl_height_bulk_richardson_pt25'][:])
        pbl_reg0.append(fi.variables['pbl_regime_type_liu_liang'][:])

        pbl_r50.append(fi.variables['pbl_height_bulk_richardson_pt5'][:])
        qc_pbl_r50.append(fi.variables['qc_pbl_height_bulk_richardson_pt5'][:])

        pbl_hef0.append(fi.variables['pbl_height_heffter'][:])
        qc_pbl_hef0.append(fi.variables['qc_pbl_height_heffter'][:])
    
        pres0.append(fi.variables['atm_pres'][:])
        temp0.append(fi.variables['air_temp'][:])
        
        pbl_time0.append(datetime(int(i[26:30]),int(i[30:32]),int(i[32:34]),int(i[35:37]),int(i[37:39]),int(i[39:41])))
        fi.close(); #nc2.close()

#                  flag_meanings: neutral_boundary_layer stable_boundary_layer convective_boundary_layer
#                  flag_0_description: Neutral boundary layer
#                  flag_1_description: Stable boundary layer
#                  flag_2_description: Convective boundary layer

# Convert to np.arrays
pres = np.asarray(list(itertools.chain(*pres0)));
temp = np.asarray(list(itertools.chain(*temp0)));
pbl_ll = np.asarray(list(itertools.chain(pbl_ll0)));
qc_pbl_ll = np.asarray(list(itertools.chain(qc_pbl_ll0)));
pbl_reg = np.asarray(list(itertools.chain(pbl_reg0)));
pbl_r25 = np.asarray(list(itertools.chain(pbl_r250)));
qc_pbl_r25 = np.asarray(list(itertools.chain(qc_pbl_r250)));
pbl_time = copy.deepcopy(pbl_time0)
pbl_r5 = np.asarray(list(itertools.chain(pbl_r50)));
qc_pbl_r5 = np.asarray(list(itertools.chain(qc_pbl_r50)));
pbl_hef = np.asarray(list(itertools.chain(pbl_hef0)));
qc_pbl_hef = np.asarray(list(itertools.chain(qc_pbl_hef0)));

# Use qc checks in files to eliminate questionable data
for i in np.arange(0,len(qc_pbl_ll)):
    if qc_pbl_ll[i] > 0:
        pbl_ll[i] = np.nan
        pbl_reg[i] = np.nan

for i in np.arange(0,len(qc_pbl_r25)):
    if qc_pbl_r25[i] > 0:
        pbl_r25[i] = np.nan
        
for i in np.arange(0,len(qc_pbl_r5)):
    if qc_pbl_r5[i] > 0:
        pbl_r5[i] = np.nan

for i in np.arange(0,len(qc_pbl_hef)):
    if qc_pbl_hef[i] > 0:
        pbl_hef[i] = np.nan

# Convert to a pandas dataframe
d = {'time': pbl_time, 'pbl_ll': pbl_ll, 'qc_pbl_ll': qc_pbl_ll, 'pbl_reg': pbl_reg, 'pbl_r25': pbl_r25, 'qc_pbl_r25': qc_pbl_r25, 'pbl_r5': pbl_r5, 'qc_pbl_r5': qc_pbl_r5, 'pbl_hef': pbl_hef, 'qc_pbl_hef': qc_pbl_hef}
df = pd.DataFrame(data=d)
df = df.set_index('time')
hhh = df.index.hour; mmm = df.index.month; yyy = df.index.year;
df['yr'] = yyy;    df['mo'] = mmm;    df['hr'] = hhh

df_3H = df.resample('180Min',loffset='90Min').mean()
new_date_index_3H = [datetime(2009,1,1,1,30,0) + timedelta(hours=3*x) for x in range(0,14600)]
df_3H = df_3H.reindex(new_date_index_3H, fill_value = np.nan)
df_3H['yr'] = df_3H.index.year; df_3H['mo'] = df_3H.index.month; df_3H['hr'] = df_3H.index.hour

df_2H = df.resample('120Min',loffset='60Min').mean()
new_date_index_2H = [datetime(2009,1,1,1,0,0) + timedelta(hours=2*x) for x in range(0,21900)]
df_2H = df_2H.reindex(new_date_index_2H, fill_value = np.nan)
df_2H['yr'] = df_2H.index.year; df_2H['mo'] = df_2H.index.month; df_2H['hr'] = df_2H.index.hour

# Group data by season
df_cur = copy.deepcopy(df_3H); nt = 8; dt = 3; mo_num = 8;
df_cur_MAM = df_cur[(df_cur.iloc[:,mo_num] == 3) | (df_cur.iloc[:,mo_num] == 4) | (df_cur.iloc[:,mo_num] == 5)]
df_cur_JJA = df_cur[(df_cur.iloc[:,mo_num] == 6) | (df_cur.iloc[:,mo_num] == 7) | (df_cur.iloc[:,mo_num] == 8)]
df_cur_SON = df_cur[(df_cur.iloc[:,mo_num] == 9) | (df_cur.iloc[:,mo_num] == 10) | (df_cur.iloc[:,mo_num] == 11)]
df_cur_D = df_cur[(df_cur.iloc[:,mo_num] == 12)]
df_cur_JF = df_cur[(df_cur.iloc[:,mo_num] == 1) | (df_cur.iloc[:,mo_num] == 2)]
df_cur_DJF = df_cur[(df_cur.iloc[:,mo_num] == 1) | (df_cur.iloc[:,mo_num] == 2) | (df_cur.iloc[:,mo_num] == 12)]

pbl_plot = np.reshape(df_cur_JJA.pbl_ll,[-1, nt])

x_arr = np.arange(dt/2,24,dt)+1
fig,ax = plt.subplots(1,1,figsize = (8,5))
ax.scatter(x_arr,np.nanmedian(pbl_plot,axis=0))
ax.scatter(x_arr,np.nanpercentile(pbl_plot,75,axis=0))
ax.scatter(x_arr,np.nanpercentile(pbl_plot,25,axis=0))
ax.set_xticks(x_arr)
ax.set_xticklabels(x_arr)
ax.set_xlabel('UTC time')
ax.set_ylabel('Planetary Boundary Layer Height (m)')

ax2 = ax.twiny()
new_lab = x_arr-5
for i in np.arange(0,len(new_lab)):
    if new_lab[i] < 0:
        new_lab[i] = new_lab[1]+24
new_pos = x_arr
ax2.set_xticks(new_pos)
ax2.set_xticklabels(new_lab)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('outward',45))
ax2.set_xlabel('Central Daylight Time')
ax2.set_xlim(ax.get_xlim())

plt.tight_layout()

for i in np.arange(0,nt):
    print(np.count_nonzero(~np.isnan(pbl_plot[:,i])))

powerpath = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/H/Spectral/'
pblmeth = ['LL','R25','R5','HEF']
pblmeth = ['R25']
pctiles = [23,77]
for j in np.arange(0,4):
    colors=['blue','crimson','goldenrod','0.75']
    cnts = 0
    dx = [-1.5, -0.5, 0.5, 1.5]
    fig,ax = plt.subplots(1,1,figsize = (6,5))
    cnt_tp = np.zeros(4)
    cnt_tnp = np.zeros(4)
    
    for season in ['MAM','JJA','SON','DJF']:
        # import pickle
        import pickle
        filename = powerpath+'times_power_'+season+'_intN0_7day.p'
        [tp,tp_data,tp_4ddata,tp_pwr,tnp,tnp_data,tnp_4ddata,tnp_pwr] = pickle.load(open(filename,"rb"))
        savename_comp = season+'_Comparison_Power_H.png'
        savename_comp2 = season+'_Comparison_Power_H.eps'
       
        cdata = copy.deepcopy(tp); cdata2 = copy.deepcopy(tp_data)
        p_num =len(cdata)
        p_reshape_vals = np.zeros((1,8))
        for i in np.arange(0,len(cdata)):
            temp_date_arr = [np.datetime64(df_time)-np.datetime64(cdata[i][0]) for df_time in df_3H.index]
            id1 = np.argmin(np.abs(temp_date_arr))
        
            temp_date_arr2 = [np.datetime64(df_time)-np.datetime64(cdata[i][-1]) for df_time in df_3H.index]
            id2 = np.argmin(np.abs(temp_date_arr2))
        
    #        savename = 'Power'+str(np.datetime64(cdata[i][0]))[0:13]  
    #        fig,ax = plt.subplots(1,1,figsize=[10,7])
    #        ax.plot(cdata[i],cdata2[i])    
    #        ax.set_ylabel('Aerosol Conc.')
    #        
    #        ax2 = ax.twinx()       
    #        ax2.scatter(df_3H.pbl_ll.iloc[id1:id2+1].index,df_3H.pbl_ll.iloc[id1:id2+1].values)
    #        ax2.scatter(df_3H.pbl_r25.iloc[id1:id2+1].index,df_3H.pbl_r25.iloc[id1:id2+1].values,color='r')
    #        ax2.scatter(df_3H.pbl_reg.iloc[id1:id2+1].index,df_3H.pbl_reg.iloc[id1:id2+1].values*100,color='y')
    #
    ##        ax2 = ax.twinx()       
    ##        ax2.scatter(df_3H.pbl_ll.iloc[id1+12:id2-11].index,df_3H.pbl_ll.iloc[id1+12:id2-11].values)
    ##        ax2.scatter(df_3H.pbl_r25.iloc[id1+12:id2-11].index,df_3H.pbl_r25.iloc[id1+12:id2-11].values,color='r')
    ##        ax2.scatter(df_3H.pbl_reg.iloc[id1+12:id2-11].index,df_3H.pbl_reg.iloc[id1+12:id2-11].values*100,color='y')
    ## 
    #        ax2.set_ylabel('Boundary Layer Height (m)')
    #        plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./PBL/'+savename+'.png')
    #        plt.close(fig)
    
            if tp_pwr[i] >= np.nanpercentile(tp_pwr+tnp_pwr,pctiles[1]):
                cnt_tp[cnts] = cnt_tp[cnts]+1
                if str(cdata[i][0])[11:13] == '01':
                    if j == 0:
                        stack_vals = np.reshape(df_3H.pbl_ll.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    elif j == 1:
                        stack_vals = np.reshape(df_3H.pbl_r25.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    elif j == 2:
                        stack_vals = np.reshape(df_3H.pbl_r5.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    elif j == 3:
                        stack_vals = np.reshape(df_3H.pbl_hef.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    stack_vals_time = np.reshape(df_3H.index[id1:id2+1].values,[7, -1])[1:-1]
                    cdata_plot = np.reshape(cdata[i],[7, -1])[1:-1]
                    cdata2_plot = np.reshape(cdata2[i],[7, -1])[1:-1]
                elif str(cdata[i][0])[11:13] == '13' :
                    if j == 0:
                        stack_vals = np.reshape(df_3H.pbl_ll.iloc[id1+12:id2-11].values,[4, -1])[:]
                    elif j == 1:
                        stack_vals = np.reshape(df_3H.pbl_r25.iloc[id1+12:id2-11].values,[4, -1])[:]
                    elif j == 2:
                        stack_vals = np.reshape(df_3H.pbl_r5.iloc[id1+12:id2-11].values,[4, -1])[:]
                    elif j == 3:
                        stack_vals = np.reshape(df_3H.pbl_hef.iloc[id1+12:id2-11].values,[4, -1])[:]
                        
                    stack_vals_time = np.reshape(df_3H.index[id1+12:id2-11].values,[4, -1])[:]            
                    cdata_plot = np.reshape(cdata[i][18:66],[4, -1])[:]
                    cdata2_plot = np.reshape(cdata2[i][18:66],[4, -1])[:]
                    
                nax = np.shape(stack_vals)[0]
                savename = pblmeth[j]+'_PowerH_'+str(tp_pwr[i])+str(np.datetime64(cdata[i][0]))[0:13]   
                fig,axs = plt.subplots(nax,1,figsize=[7,int(nax*2)])
                for iii in np.arange(0,nax):
                    axs[iii].plot(cdata_plot[iii,:],cdata2_plot[iii,:])    
                    axs[iii].set_ylabel('Aerosol Conc.')
                
                    ax2 = axs[iii].twinx()       
                    ax2.scatter(stack_vals_time[iii,:],stack_vals[iii,:],color='b')
                    ax2.set_ylabel('Boundary Layer Height (m)')
                axs[0].set_title(str(tp_pwr[i]))
                plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./PBL/x'+savename+'.png')
                plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./PBL/x'+savename+'.eps')
                plt.close(fig)
                    
                p_reshape_vals = np.vstack((p_reshape_vals,stack_vals))
        
                        
        cdata = copy.deepcopy(tnp); cdata2 = copy.deepcopy(tnp_data)
        np_reshape_vals = np.zeros((1,8))
        np_num =len(cdata)
        for i in np.arange(0,len(cdata)):
            temp_date_arr = [np.datetime64(df_time)-np.datetime64(cdata[i][0]) for df_time in df_3H.index]
            id1 = np.argmin(np.abs(temp_date_arr))
        
            temp_date_arr2 = [np.datetime64(df_time)-np.datetime64(cdata[i][-1]) for df_time in df_3H.index]
            id2 = np.argmin(np.abs(temp_date_arr2))
        
            if tnp_pwr[i] <= np.nanpercentile(tp_pwr+tnp_pwr,pctiles[0]):
                cnt_tnp[cnts] = cnt_tnp[cnts]+1
                if str(cdata[i][0])[11:13] == '01':
                    if j == 0:
                        stack_vals = np.reshape(df_3H.pbl_ll.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    elif j == 1:
                        stack_vals = np.reshape(df_3H.pbl_r25.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    elif j == 2:
                        stack_vals = np.reshape(df_3H.pbl_r5.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    elif j == 3:
                        stack_vals = np.reshape(df_3H.pbl_hef.iloc[id1:id2+1].values,[7, -1])[1:-1]
                    stack_vals_time = np.reshape(df_3H.index[id1:id2+1].values,[7, -1])[1:-1]
                    cdata_plot = np.reshape(cdata[i],[7, -1])[1:-1]
                    cdata2_plot = np.reshape(cdata2[i],[7, -1])[1:-1]
                elif str(cdata[i][0])[11:13] == '13' :
                    if j == 0:
                        stack_vals = np.reshape(df_3H.pbl_ll.iloc[id1+12:id2-11].values,[4, -1])[:]
                    elif j == 1:
                        stack_vals = np.reshape(df_3H.pbl_r25.iloc[id1+12:id2-11].values,[4, -1])[:]
                    elif j == 2:
                        stack_vals = np.reshape(df_3H.pbl_r5.iloc[id1+12:id2-11].values,[4, -1])[:]
                    elif j == 3:
                        stack_vals = np.reshape(df_3H.pbl_hef.iloc[id1+12:id2-11].values,[4, -1])[:]

                    stack_vals_time = np.reshape(df_3H.index[id1+12:id2-11].values,[4, -1])[:]            
                    cdata_plot = np.reshape(cdata[i][18:66],[4, -1])[:]
                    cdata2_plot = np.reshape(cdata2[i][18:66],[4, -1])[:]

                nax = np.shape(stack_vals)[0]
                savename = pblmeth[j]+'_NoPowerH_'+str(tnp_pwr[i])+str(np.datetime64(cdata[i][0]))[0:13]   
                fig,axs = plt.subplots(nax,1,figsize=[7,int(nax*2)])
                for iii in np.arange(0,nax):
                    axs[iii].plot(cdata_plot[iii,:],cdata2_plot[iii,:])    
                    axs[iii].set_ylabel('Aerosol Conc.')
                
                    ax2 = axs[iii].twinx()       
                    ax2.scatter(stack_vals_time[iii,:],stack_vals[iii,:],color='b')
#                    ax2.scatter(df_3H.pbl_reg.iloc[id1:id2+1].index,df_3H.pbl_reg.iloc[id1:id2+1].values*100,color='y')
                    ax2.set_ylabel('Boundary Layer Height (m)')
                axs[0].set_title(str(tnp_pwr[i]))
                plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./PBL/x'+savename+'.png')
                plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./PBL/x'+savename+'.eps')
                plt.close(fig)
    ##        
                np_reshape_vals = np.vstack((np_reshape_vals,stack_vals))
        
        ids = [1,3,5,7]
        ax.scatter(x_arr[ids]+0.25+dx[cnts],np.nanpercentile(p_reshape_vals,25,axis=0)[ids],marker='_',color=colors[cnts],label='_nolegend_')
        ax.scatter(x_arr[ids]+0.25+dx[cnts],np.nanpercentile(p_reshape_vals,50,axis=0)[ids],s=50,marker='o',color=colors[cnts],edgecolors='w',label=season+' High Power ('+str(int(cnt_tp[cnts]))+')')
        ax.scatter(x_arr[ids]+0.25+dx[cnts],np.nanpercentile(p_reshape_vals,75,axis=0)[ids],marker='_',color=colors[cnts],label='_nolegend_')
    
        ax.scatter(x_arr[ids]-0.25+dx[cnts],np.nanpercentile(np_reshape_vals,25,axis=0)[ids],marker='_',color=colors[cnts],label='_nolegend_')
        ax.scatter(x_arr[ids]-0.25+dx[cnts],np.nanpercentile(np_reshape_vals,50,axis=0)[ids],s=50,marker='d',color=colors[cnts],edgecolors='w',label=season+' Low Power ('+str(int(cnt_tnp[cnts]))+')')
        ax.scatter(x_arr[ids]-0.25+dx[cnts],np.nanpercentile(np_reshape_vals,75,axis=0)[ids],marker='_',color=colors[cnts],label='_nolegend_')
        
    
    
    #    exec('pbl_plot=np.reshape(df_cur_'+season+'.pbl_ll,[-1, nt])')
    #    exec('pbl_plot=np.reshape(df_cur_'+season+'.pbl_r25,[-1, nt])')
    #    ax.scatter(x_arr[ids],np.nanpercentile(pbl_plot,25,axis=0)[ids],marker='^',color='gray',label='_nolegend_')
    #    ax.scatter(x_arr[ids],np.nanpercentile(pbl_plot,50,axis=0)[ids],marker='o',color='gray',label='All Season')
    #    ax.scatter(x_arr[ids],np.nanpercentile(pbl_plot,75,axis=0)[ids],marker='v',color='gray',label='_nolegend_')
    
        ax.set_ylim([0,3000])
        ax.set_xticks(x_arr[ids])
        ax.set_xticklabels(['05:30','11:30','17:30','23:30'])
        ax.set_xlabel('UTC Time')
        ax.set_ylabel('Boundary Layer Height (m)')
        #ax.legend()
        
        
        cnts = cnts + 1
    
    ax.set_ylim([0,3000])
    ax2 = ax.twiny()
    new_lab = x_arr-5
    for i in np.arange(0,len(new_lab)):
        if new_lab[i] < 0:
            new_lab[i] = new_lab[1]+24
    new_pos = x_arr
    ax2.set_xticks(new_pos[ids])
    ax2.set_xticklabels(['00:30','06:30','12:30','18:30'])
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.spines['bottom'].set_position(('outward',45))
    ax2.set_xlabel('Central Daylight Time')
    ax2.set_xlim(ax.get_xlim())

    ax.legend(loc=2,fontsize=10)
    
    plt.tight_layout()
    plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/PBL/'+pblmeth[j]+'_'+str(pctiles[0])+str(pctiles[1])+'_equal_wIQR_1FIG_'+savename_comp)
    plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/PBL/'+pblmeth[j]+'_'+str(pctiles[0])+str(pctiles[1])+'_equal_wIQR_1FIG_'+savename_comp2)
    plt.close(fig)

    
pathname = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/'
[df_120M,DP_full2] = pickle.load(open(pathname+'df_120M_FnvFinal_D2E4_var95_LT18_test2.p',"rb"))

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

N0_DJF_res = np.reshape(df_120M_DJF.intN0,[-1, 12])
N0_times = np.reshape(df_120M_DJF.intN0,[-1, 12])

weekly_intN0 = df_120M.intN0.values
weekly_intN0 = np.reshape(weekly_intN0,[-1, 12])
weekly_intN0[weekly_intN0 == 0] = np.nan
dates = df_120M.intN0.index.values
dates = np.reshape(dates,[-1, 12])
new_dates = np.zeros((1,1))
new_data = np.zeros((1,12))
for i in np.arange(0,np.shape(weekly_intN0)[0]):
    if np.isnan(weekly_intN0[i,:]).any() == False:
        new_dates = np.vstack((new_dates,i))
        new_data = np.vstack((new_data,weekly_intN0[i,:]))
new_data = new_data[1:,:]
new_dates = new_dates[1:]

df_dates = []
for i in np.arange(0,len(new_dates)):
    df_dates.append(dates[int(new_dates[i,0])][0])
    
max_id = []
conc_mult = []
for i in np.arange(0,len(new_dates)):
    max_id.append(np.argmax(new_data[i,:]))
    conc_mult.append(new_data[i,max_id[i]]/np.mean(new_data[i,:]))

for i in np.arange(0,len(new_dates)):
    if conc_mult[i] > 4:       
        temp_date_arr = [np.datetime64(df_time)-np.datetime64(df_dates[i]) for df_time in df_3H.index]
        id1 = np.argmin(np.abs(temp_date_arr))
        
        
        fig,ax = plt.subplots(1,1,figsize=(8,5))
        ax.plot(dates[int(new_dates[i,0])],new_data[i,:])
            
        ax2 = ax.twinx()
        ax2.scatter(df_3H.index[id1:id1+8],df_3H.pbl_ll.iloc[id1:id1+8])
        ax2.set_xlim(ax.get_xlim())
        ax.set_title(str(df_dates[i])[0:10])

        ax.set_ylabel('Aerosol Concentration')
        ax2.set_ylabel('Boundary Layer Height (m)')
        plt.tight_layout()
        plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/PBL/Aerosol_PBL_Daily_'+str(df_dates[i])[0:10]+'.png')
        plt.savefig('/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Paper/Final_Submission/Reviews/Final_Accepted/PBL/Aerosol_PBL_Daily_'+str(df_dates[i])[0:10]+'.eps')
        plt.close(fig)
