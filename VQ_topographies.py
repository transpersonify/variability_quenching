#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 00:39:46 2021

@author: sn254804
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 17:48:56 2021

@author: sn254804
"""

import scipy.io as spio
import numpy as np
from scipy.signal import hilbert
import matplotlib.pyplot as plt
import mne
import os
import scipy.stats as spstats
import pickle
import seaborn as sns
from scipy.spatial.distance import pdist
montage = mne.channels.make_standard_montage('GSN-HydroCel-256')
import multiprocessing
from functools import partial
import glob

info = mne.create_info(montage.ch_names,250,ch_types='eeg')
#info['bads']  = bads
info.set_montage(montage)
con,ch_names = mne.channels.find_ch_adjacency(info, "eeg")
neighbour_list = []

C = 256
for i in range(C):
    neighbour_list.append(np.where(con.toarray()[i])[0])
    

# =============================================================================
# ages = spio.loadmat('/volatile/home/sn254804/Bureau/script_panel_generator/ages/ages.mat')
# ages = np.squeeze(ages['age'])
# 
# young = np.where(ages<=12)[0]
# old = np.where(ages>=16)[0]
# 
# =============================================================================

def cov_topo(wd,files,neighbour_list,sub): #, keepChannels)
    data = spio.loadmat(files[sub])
    #keepChannels = np.squeeze(keepChn[sub]) - 1
    T = 475
    covL_all = np.zeros((C,T))
    covR_all = np.zeros((C,T))
    print('Now running... ' + files[sub])
    for idn,neigh in enumerate(neighbour_list):
        
        print("Evaluating Cluster..." + str(idn))
        #erp_clusters = np.intersect1d(keepChannels,neigh)
        erp_clusters = neigh
        left_trials = data['left_trials'][erp_clusters,:,:]
        right_trials = data['right_trials'][erp_clusters,:,:]
# =============================================================================
#         C = len(erp_clusters)
#     
#         Nl = left_trials.shape[2]
#         Nr = right_trials.shape[2]
#         
# =============================================================================
        T = left_trials.shape[1]
            
        for t in range(T):
            covL_all[idn,t] = np.nanmean(pdist(left_trials[:,t,:].T,'correlation'))
            covR_all[idn,t] = np.nanmean(pdist(right_trials[:,t,:].T,'correlation'))
       
    return covL_all,covR_all


if __name__ == "__main__":
    #channels = spio.loadmat('/volatile/home/sn254804/Bureau/script_panel_generator/Bad_channels/lateral_faces/keepChannels.mat')
    #keepChn = np.squeeze(channels['keepChn'])
    N = 13
    p = multiprocessing.Pool(6)
    wd = '//volatile/home/sn254804/Bureau/phase_reset_analysis/DATA/ADULTS/'
    files = glob.glob(wd + '*.mat')
    files.sort()
    files = files[:-1]
    func = partial(cov_topo,wd,files,neighbour_list)
    #[covL_all,covR_all] = cov_topo(wd,files,keepChn,neighbour_list,2)
    data_cov = p.map(func, np.arange(N))

path = '/volatile/home/sn254804/Bureau/phase_reset_analysis/DATA/ADULTS/VQ_topography/'



for i in range(13):
    [covL,covR] = data_cov[i]
    data_dict = {'covL_topo': covL, 'covR_topo': covR}
    spio.savemat(path + 'VQ_topo_'+ files[i].split('/')[-1],data_dict)


for i in range(13):
    [covL,covR] = data_cov[i]
    data_dict = {'covL_topo': covL, 'covR_topo': covR}
    spio.savemat(path + 'VQ_topo_'+ files[i].split('/')[-1],data_dict)

# VQ_files = '/volatile/home/sn254804/Bureau/script_panel_generator/data/VQ/Infants/lateral_faces/'
# with open(VQ_files + 'VQ_no_filter.pkl','rb') as ff:
#     [left,right] = pickle.load(ff)
# leftvq = spstats.zscore(left,axis=1)
# rightvq = spstats.zscore(right,axis=1)
# times = np.arange(-0.4,1.5,0.004)

# CV_files = '/volatile/home/sn254804/Bureau/phase_reset_analysis/'
# with open(CV_files + 'Infants_CV_all.pkl','rb') as ff:
#     [leftCV,rightCV] = pickle.load(ff)

