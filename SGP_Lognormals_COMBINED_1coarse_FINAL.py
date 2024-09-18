#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 22:22:26 2018

@author: petermarinescu
"""

#### Need to:
import pickle
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import numpy as np

filename = '/Volumes/GoogleDrive/My Drive/Research_Projects/CCN_Closure_M.S./Size_Distributions_Sam_Reorder.p'
#filename = '/Volumes/Data/Act_Spectrum/Final/Size_Distributions_Sam_Reorder.p'
#filename = '/Volumes/Data/Act_Spectrum/H/Median_Size_Distributions.p'

test = pickle.load(open(filename,"rb"))

def make_lognormal(xvals,d_med,std_g,N_t):
    import numpy as np
    i=0; yvals = np.zeros(len(xvals))
    for d in xvals:
        a = (N_t / (np.log(std_g)*np.sqrt(2*np.pi)));
        b = -(np.power(np.log(d)-np.log(d_med),2.0))/(2*np.power(np.log(std_g),2.0));
        yvals[i] = a * np.exp(b); 
        i = i + 1;
    return yvals
        
diams = test[0]
diam_mids = (diams[0:len(diams)-1]+diams[1:])/2
diam_mids = np.concatenate((diam_mids[0]-(np.diff(diam_mids[0:2])),diam_mids,diam_mids[-1]+(np.diff(diam_mids[-2:]))),axis=0)
logDP_mids = np.log10(diam_mids)
dlogDP_mids = np.diff(logDP_mids)

lnDP_mids = np.log(diam_mids)
dlnDP_mids = np.diff(lnDP_mids)


seasons = ['ALL','MAM','JJA','SON','DJF']

v_scale0 = 1.05
v_scale1 = 1.075
v_scale2 = 1.10
# Initialize Arrays
r_n = np.zeros((5,5)); d_n = np.zeros((5,5)); s_n = np.zeros((5,5)); t_n = np.zeros((5,5))
r_s = np.zeros((5,5)); d_s = np.zeros((5,5)); s_s = np.zeros((5,5)); t_s = np.zeros((5,5))
r_v = np.zeros((5,5)); d_v = np.zeros((5,5)); s_v = np.zeros((5,5)); t_v = np.zeros((5,5))

# MAM Volume modes from Sam
# Volume n smallest number modes from Sam
d_v[1,1] = 0.156; s_v[1,1] = 1.81; t_v[1,1] = 1.32 * np.log10(np.exp(1))
d_v[2,1] = 0.279; s_v[2,1] = 1.56; t_v[2,1] = 4.25  * np.log10(np.exp(1))

d_v[3,1] = 2.15; s_v[3,1] = 1.68; t_v[3,1] = 2.9  * np.log10(np.exp(1))
d_v[4,1] = 5.95; s_v[4,1] = 1.44; t_v[4,1] = 1.1  * np.log10(np.exp(1))

d_v[3,1] = 2.88; s_v[3,1] = 1.99; t_v[3,1] = 4.05*v_scale0  * np.log10(np.exp(1))
d_v[4,1] = 0; s_v[4,1] = 1; t_v[4,1] = 0.0

d_n[0,1] = 0.0055; s_n[0,1] = 2.8; t_n[0,1] = 7100  * np.log10(np.exp(1))
r_n[0,1] = d_n[0,1] / 2;

#d_n[0,1] = 0.0045; s_n[0,1] = 2.8; t_n[0,1] = 7500
#r_n[0,1] = d_n[0,1] / 2;

# JJA Volume modes from Sam
# Volume n smallest number modes from Sam
d_v[1,2] = 0.1685; s_v[1,2] = 1.76; t_v[1,2] = 2.805  * np.log10(np.exp(1))
d_v[2,2] = 0.293; s_v[2,2] = 1.56; t_v[2,2] = 5.63  * np.log10(np.exp(1))

d_v[3,2] = 2.47; s_v[3,2] = 1.60; t_v[3,2] = 4.3 * np.log10(np.exp(1))
d_v[4,2] = 6.; s_v[4,2] = 1.4; t_v[4,2] = 1.5 * np.log10(np.exp(1))

d_v[3,2] = np.exp(((t_v[4,2]/(t_v[3,2]+t_v[4,2]))*np.log(d_v[4,2])+(t_v[3,2]/(t_v[3,2]+t_v[4,2]))*np.log(d_v[3,2])))*1.01; 
s_v[3,2] = 1.93; t_v[3,2] = (t_v[3,2]+t_v[4,2])*(v_scale0-0.04)
d_v[4,2] = 0; s_v[4,2] = 1; t_v[4,2] = 0

d_n[0,2] = 0.0055; s_n[0,2] = 2.8; t_n[0,2] = 5000  * np.log10(np.exp(1))
#d_n[0,2] = 0.0045; s_n[0,2] = 2.8; t_n[0,2] = 5400  * np.log10(np.exp(1))

r_n[0,2] = d_n[0,2] / 2;

# SON Volume modes from Sam
# Volume n smallest number modes from Sam
d_v[1,3] = 0.148; s_v[1,3] = 1.78; t_v[1,3] = 1.66  * np.log10(np.exp(1))
d_v[2,3] = 0.273; s_v[2,3] = 1.54; t_v[2,3] = 4.15  * np.log10(np.exp(1))

d_v[3,3] = 2.176; s_v[3,3] = 1.58; t_v[3,3] = 1.85  * np.log10(np.exp(1))
d_v[4,3] = 6.041; s_v[4,3] = 1.56; t_v[4,3] = 2.0  * np.log10(np.exp(1))

d_v[3,3] = np.exp(((t_v[4,3]/(t_v[3,3]+t_v[4,3]))*np.log(d_v[4,3])+(t_v[3,3]/(t_v[3,3]+t_v[4,3]))*np.log(d_v[3,3]))); 
s_v[3,3] = 2.0; t_v[3,3] = (t_v[3,3]+t_v[4,3])*v_scale0  
d_v[4,3] = 0; s_v[4,3] = 1; t_v[4,3] = 0

d_n[0,3] = 0.0055; s_n[0,3] = 2.8; t_n[0,3] = 6700  * np.log10(np.exp(1))
r_n[0,3] = d_n[0,3] / 2;

# DJF Volume modes from Sam
# Volume n smallest number modes from Sam
d_v[1,4] = 0.163; s_v[1,4] = 1.84; t_v[1,4] = 1.89  * np.log10(np.exp(1))
d_v[2,4] = 0.302; s_v[2,4] = 1.54; t_v[2,4] = 5.2  * np.log10(np.exp(1))

d_v[3,4] = 1.94; s_v[3,4] = 1.56; t_v[3,4] = 1.36  * np.log10(np.exp(1))
d_v[4,4] = 5.276; s_v[4,4] = 1.45; t_v[4,4] = 1.29  * np.log10(np.exp(1))

d_v[3,4] = np.exp(((t_v[3,4]/(t_v[3,4]+t_v[4,4]))*np.log(d_v[4,4])+(t_v[3,4]/(t_v[3,4]+t_v[4,4]))*np.log(d_v[3,4]))); 
s_v[3,4] = 1.94; t_v[3,4] = (t_v[3,4]+t_v[4,4])*v_scale2
d_v[4,4] = 0; s_v[4,4] = 0; t_v[4,4] = 0

d_n[0,4] = 0.0045; s_n[0,4] = 2.8; t_n[0,4] = 4400  * np.log10(np.exp(1))
r_n[0,4] = d_n[0,4] / 2;

# ALL Volume modes from Sam
# Volume n smallest number modes from Sam
d_v[1,0] = 0.172; s_v[1,0] = 1.82; t_v[1,0] = 2.3  * np.log10(np.exp(1))
d_v[2,0] = 0.286; s_v[2,0] = 1.53; t_v[2,0] = 4.4  * np.log10(np.exp(1))

d_v[3,0] = 2.4; s_v[3,0] = 1.7; t_v[3,0] = 2.83  * np.log10(np.exp(1))
d_v[4,0] = 6.148; s_v[4,0] = 1.43; t_v[4,0] = (1.23*1.05)*v_scale1 * np.log10(np.exp(1))

d_v[3,0] = np.exp(((t_v[4,0]/(t_v[3,0]+t_v[4,0]))*np.log(d_v[4,0])+(t_v[3,0]/(t_v[3,0]+t_v[4,0]))*np.log(d_v[3,0]))); 
s_v[3,0] = 1.97; t_v[3,0] = t_v[3,0]+t_v[4,0]
d_v[4,0] = 0; s_v[4,0] = 1; t_v[4,0] = 0

d_n[0,0] = 0.0053; s_n[0,0] = 2.8; t_n[0,0] = 6000  * np.log10(np.exp(1))
r_n[0,0] = d_n[0,0] / 2;

modes_v = np.zeros((len(diams),5,len(seasons)))
modes_s = np.zeros((len(diams),5,len(seasons)))
modes_n = np.zeros((len(diams),5,len(seasons)))

for j in np.arange(0,len(seasons)):
    
    # Convert diameters to radii
    for i in np.arange(0,5):
        r_v[i,j] = d_v[i,j] / 2

    # Convert modes 1-4 volume data to number data
    for i in np.arange(1,5):    
        r_n[i,j] = np.exp(np.log(r_v[i,j])-(3*np.power(np.log(s_v[i,j]),2.0)))
        d_n[i,j] = r_n[i,j]*2; print(d_v[i,j])
        s_n[i,j] = s_v[i,j]
        t_n[i,j] = (t_v[i,j]*3/np.pi/4) / np.exp(3*np.log(r_n[i,j])+4.5*(np.power(np.log(s_v[i,j]),2)))

    # convert number data to volume and surface area
    for i in np.arange(0,5):
        r_v[i,j] = np.exp(np.log(r_n[i,j])+(3*np.power(np.log(s_n[i,j]),2.0)))
        s_v[i,j] = s_n[i,j]
        t_v[i,j] = t_n[i,j]*4*np.pi/3*np.exp(3*np.log(r_n[i,j])+4.5*np.power(np.log(s_n[i,j]),2.0))
        d_v[i,j] = r_v[i,j]*2
        
        r_s[i,j] = np.exp(np.log(r_n[i,j])+(2*np.power(np.log(s_n[i,j]),2.0)))
        s_s[i,j] = s_n[i,j]
        t_s[i,j] = t_n[i,j]*4*np.pi*np.exp(2*np.log(r_n[i,j])+2*np.power(np.log(s_n[i,j]),2.0))
        d_s[i,j] = r_s[i,j]*2

    # Make Modes   
    for i in np.arange(0,4):
        modes_v[:,i,j] = make_lognormal(diams,r_v[i,j]*2,s_v[i,j],t_v[i,j])
        modes_s[:,i,j] = make_lognormal(diams,r_s[i,j]*2,s_s[i,j],t_s[i,j])
        modes_n[:,i,j] = make_lognormal(diams,r_n[i,j]*2,s_n[i,j],t_n[i,j])

#for i in np.arange(0,5):
#    for j in np.arange(0,5):
#        print(r_n[i,j])

#for i in np.arange(0,5):
#    for j in np.arange(0,5):
#        print(str(np.nansum(modes_n[:,i,j]*dlnDP_mids))+' | '+str(t_n[i,j])+' | '+str(d_n[i,j]))

colors=['0.15','0.35','0.55','0.70']
lw_m = 4
lw_i = 2

N_dif = np.zeros(5); N_pdif = np.zeros(5);
S_dif = np.zeros(5); S_pdif = np.zeros(5);
V_dif = np.zeros(5); V_pdif = np.zeros(5);

fs = 16.5
plt.rcParams.update({'font.size':15})
fig,ax = plt.subplots(3,5,figsize=(17.9,7.5),sharex=True)

bin_edges = [0, 0.029, 0.140, 0.800, 30] # bin edges in microns
max_vals = [2500,120,5]
for j in np.arange(0,len(seasons)):
    for k in np.arange(0,3):
        for m in np.arange(1,5):
            ax[k,j].plot([bin_edges[m],bin_edges[m]],[0,max_vals[k]],'-',lw = 2.5,c='0.5')

for j in np.arange(0,len(seasons)):

    sum_mode_n = np.zeros(len(diams)); sum_mode_s = np.zeros(len(diams)); sum_mode_v = np.zeros(len(diams))
    # Plot Modes    
    ax[0,j].plot(diams,test[1][:,j],'-k',lw=lw_m,label='SGP Data')
    for i in np.arange(0,4):
        sum_mode_n[:] = sum_mode_n[:]+modes_n[:,i,j] 
    ax[0,j].plot(diams,sum_mode_n,color = 'goldenrod',ls=':',lw=lw_m,label='Modes Combined')
    for i in np.arange(0,4):
        ax[0,j].plot(diams,modes_n[:,i,j],color=colors[i],ls='--',lw=lw_i,label='Mode '+str(i+1))
#    ax[0,j].set_xscale('log')
    ax[0,0].set_ylabel('dN dlnD$_{p}^{-1}$ \n (# cm$^{-3}$)',fontsize=fs)
    ax[0,j].set_title(seasons[j])
    ax[0,j].set_ylim([0,2000])
    ax[0,j].set_xlim([0.007,14])

    ax[1,j].plot(diams,test[2][:,j],'-k',lw=lw_m,label='SGP Data')
    for i in np.arange(0,4):
        sum_mode_s[:] = sum_mode_s[:]+modes_s[:,i,j] 
    ax[1,j].plot(diams,sum_mode_s,color = 'goldenrod',ls=':',lw=lw_m,label='Modes Combined')
    for i in np.arange(0,4):
        ax[1,j].plot(diams,modes_s[:,i,j],color=colors[i],ls='--',lw=lw_i,label='Mode '+str(i+1))
#    ax[1,j].set_xscale('log')
    ax[1,0].set_ylabel('dSA dlnD$_{p}^{-1}$ \n ($\mu$m$^{2}$ cm$^{-3}$)',fontsize=fs)
    ax[1,j].set_ylim([0,80])
    ax[1,j].set_xlim([0.007,14])

    ax[2,j].plot(diams,test[3][:,j],'-k',lw=lw_m,label='SGP Data')
    for i in np.arange(0,4):
        sum_mode_v[:] = sum_mode_v[:]+modes_v[:,i,j] 
    ax[2,j].plot(diams,sum_mode_v,color = 'goldenrod',ls=':',lw=lw_m,label='Modes Combined')
    for i in np.arange(0,4):
        ax[2,j].plot(diams,modes_v[:,i,j],color=colors[i],ls='--',lw=lw_i,label='Mode '+str(i+1))
    ax[2,0].set_ylabel('dV dlnD$_{p}^{-1}$ \n ($\mu$m$^{3}$ cm$^{-3}$)',fontsize=fs)
    ax[2,j].set_xlabel('Diameter ($\mu$m)',fontsize=fs)
    ax[2,j].set_ylim([0,4])
    ax[2,j].set_xlim([0.007,14])
    ax[2,j].set_xscale('log')


#    ax[2,j].set_xlim([0.007,0.6])

    if j == 2:
        lgd = ax[2,j].legend(loc='lower center',bbox_to_anchor=(0.5,-0.8),ncol=7,frameon=False,fontsize=17)


    intN_orig = np.nansum(test[1][:,j]*dlnDP_mids)
    intN_fit = np.nansum(sum_mode_n*dlnDP_mids)

    N_dif[j] = intN_fit-intN_orig
    print(intN_fit)
    N_pdif[j] = N_dif[j] / intN_orig * 100

    intS_orig = np.nansum(test[2][:,j]*dlnDP_mids)
    intS_fit = np.nansum(sum_mode_s*dlnDP_mids)

    S_dif[j] = intS_fit-intS_orig
    S_pdif[j] = S_dif[j] / intS_orig * 100

    intV_orig = np.nansum(test[3][:,j]*dlnDP_mids)   
    intV_fit = np.nansum(sum_mode_v*dlnDP_mids)
    V_dif[j] = intV_fit-intV_orig
    V_pdif[j] = V_dif[j] / intV_orig * 100

        
        
#dsfs
plt.tight_layout()
fig.subplots_adjust(bottom=0.15)
plt.savefig('/Volumes/Data/Act_Spectrum/H/SGP_Distributions_Fitted_wLines_1CoarseMode_Final.png')
plt.savefig('/Volumes/Data/Act_Spectrum/H/SGP_Distributions_Fitted_wLines_1CoarseMode_Final.eps')

print(d_n[1,:])
print(s_n[1,:])

print(d_n[:,3])