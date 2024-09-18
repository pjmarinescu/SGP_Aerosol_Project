#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:30:30 2020

@author: pmarin
"""


import numpy as np
import pickle
import matplotlib.pyplot as plt
from collections import OrderedDict

names = ['i_v2','i_v3','i_f','noi_v2','noi_v3','no_f']


lt_sf = OrderedDict()
gintN = OrderedDict()
cpcint_dif = OrderedDict()
cpcint_pdif = OrderedDict()
lt_gcpc = OrderedDict()

p = '/Users/pmarin/Desktop/'
f = ['vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072120_i_v2CCNCSA.p',
     'vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072120_i_v3CCNCSA.p',
     'vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072720_i.p',
     'vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072120_noi_v2CCNCSA.p',
     'vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072120_noi_v3CCNCSA.p',
     'vFinal_D2E4_var95_LT2wk18_CPCQC_YesCCN_medpctdif_RP_072720_noi.p']


for i in np.arange(0,6):
    [lt_sf[names[i]],gintN[names[i]],cpcint_dif[names[i]],cpcint_pdif[names[i]],lt_gcpc[names[i]]] = pickle.load(open(p+f[i],"rb"))


colors = ['g','r','y','m','c','b']
ls = ['-','--',':','-','--',':']
lw = [6,4,2,6,4,2]
xlims =[10900,11000]

fig = plt.figure(figsize=[6,3])
for i in np.arange(0,6):
    plt.plot(lt_sf[names[i]],c=colors[i],ls=ls[i],lw=lw[i],label=names[i])
plt.legend()
plt.xlim(xlims)
plt.xlabel('Points')
plt.ylabel('Long-Term Scale Factor')
plt.tight_layout()


fig = plt.figure(figsize=[6,3])
for i in np.arange(0,6):
    plt.plot(cpcint_dif[names[i]],c=colors[i],ls=ls[i],lw=lw[i],label=names[i])
plt.legend()
plt.ylabel('CPCINT_DIF')
plt.xlim(xlims)
plt.ylabel('Points')
plt.tight_layout()


