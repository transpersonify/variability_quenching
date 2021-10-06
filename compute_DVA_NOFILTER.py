#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 18:26:39 2021

@author: sn254804
"""
import os
from scipy.spatial.distance import pdist
import scipy.io as spio
import pickle
import numpy as np

datapath = '/volatile/home/sn254804/Bureau/script_panel_generator/data/DVA/'
with open(datapath+ 'ADULTS_NO_FILTERING_DVA.pkl', 'rb') as ff:
    [left_dva,right_dva] = pickle.load(ff)

distpath = '/volatile/home/sn254804/Bureau/script_panel_generator/data/flybys/Adults/'
with open(distpath + 'ADULT_dist2temp_lateral_NOFILTER.pkl','rb') as ff:
    [D_left,D_right] = pickle.load(ff)

n_clus = 2
T = 200
C = 256
N = 14
fb_left_dva = [[[] for j in range(2)] for i in range(N)]
fb_right_dva = [[[] for j in range(2)] for i in range(N)]



for sub in range(13): 
    for clus in range(2):
        
        t_lim = np.arange(88,288)
        
        dist_left = D_left[sub][clus,:,:]
        dist_right = D_right[sub][clus,:,:]
        thr = np.percentile(np.vstack([dist_left[:,t_lim],dist_right[:,t_lim]]),5)
        (trial_left,T_left) = np.where(dist_left<=thr)
        (trial_right,T_right) = np.where(dist_right<=thr)
        
        
        for nl in range(trial_left.shape[0]):
            if(T_left[nl]>=t_lim[0] and T_left[nl]<=t_lim[-1]):
                t = np.arange(T_left[nl]-100,T_left[nl]+100)
                if(t[0]>=0 and t[-1]<450):
                    fb_left_dva[sub][clus].append(left_dva[sub][trial_left[nl],t])
        
        for nr in range(trial_right.shape[0]):
            if(T_right[nr]>=t_lim[0] and T_right[nr]<=t_lim[-1]):
                t = np.arange(T_right[nr]-100,T_right[nr]+100)
                if(t[0]>=0 and t[-1]<450):
                    fb_right_dva[sub][clus].append(right_dva[sub][trial_right[nr],t])
        
with open(datapath + 'ADULTS_lateral_flyby_triggerred_DVA_NO_FILTER.pkl','wb') as ff:
    pickle.dump([fb_left_dva, fb_right_dva],ff)