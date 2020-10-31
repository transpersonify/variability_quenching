# Import required packages
import scipy.io as spio
import os
import numpy as np
from scipy.spatial.distance import cdist,pdist
import warnings

def flyby_dist(templatefile,datapath,resultpath):
    """ Calculates and saves conditionwise fly-bys i.e. correlation distance between templates and each trial for each conditions separately.

        INPUTS:
        templatefile: '.mat' file containing templates.
        each variable in '.mat' file contains one array of size of channel x 1

        datapath: full path to where the '.mat' file for each subject is stored.
        Each '.mat' file contains array of size: channel x timepoints x trial

        resultpath: full path to store the distances for each subject.

        OUTPUTS:
        D : combined fly-by distance for each subject and conditions
        len(D) = #subjects,
        Each entry contains a dictionary where key = condition, value = distances
        where, distances =  N_templates x N_trials_per_cond x time array

        keys : Names for templates

        Additionally, this function will save conditionwise distances in mat file per subject in the resultpath
    """

    # Check if templatefile exists
    if(not os.path.isfile(templatefile)):
        print('ERROR: templates file does not exist \n Check if the templatepath is correct.')
        return

    # Check if datapath exists
    if(not os.path.isdir(datapath)):
        print('ERROR: Data path does not exist...')
        return

    # Check if result path exists
    if(not os.path.isdir(resultpath)):
        print('ERROR: Result path does not exist \n Create a directory to store results.')
        return

    # Combine all templates in variable CC
    templates = spio.loadmat(templatefile)
    keyvals = []
    keys = []
    args = spio.whosmat(templatefile)
    for k in range(len(args)):
        keys.append(args[k][0])
        keyvals.append(templates[args[k][0]])

    CC = np.squeeze(np.asarray(keyvals))
    n_clus = CC.shape[0]

    # List subject files and sort based on names.
    files = os.listdir(datapath)
    files.sort()
    N = len(files)      # Number of Subjects

    D = [[] for i in range(N)]  #Fly-by Distance

    # For each subject, calculate fly-by distances
    for sub in range(N):
        name = files[sub]
        print('processing subject...' , name)
        data = spio.loadmat(os.path.join(datapath, name))        # Load matrix file

        # List all conditions
        all_vars = spio.whosmat(os.path.join(datapath, name))
        conds = [all_vars[i][0] for i in range(len(all_vars))]

        dist = {}

        # For each experimental condition, calculate fly-by distance

        for args,count in zip(conds,range(len(conds))):
            trials = data[args]                     # List all trials of condition in 'args'
            if trials.size == 0:
               print('Condition: ' + conds[count] + ' is empty, keeping distance array empty...\n'),
               dist[conds[count]] = []

            else:
                if trials.ndim ==2:
                    print('Condition:' + conds[count]+ ' has only one trial, adding an extra dimension...\n'),
                    trials = trials[:,:,None]
                gfp = np.std(trials,0)               # Calculate GFP
                trials = np.divide(trials,gfp)      # Divide by GFP

                N_trials = trials.shape[2]
                T = trials.shape[1]

                d = np.zeros((n_clus,N_trials,T))
                for t in range(T):
                    d[:,:,t] = cdist(trials[:,t,:].T,CC,metric='correlation').T
                dist[conds[count]] = d
        D[sub] = dist
        spio.savemat(os.path.join(resultpath, name[:-4] + '_flybys.mat'),dist)
    return D,keys

def between_trial_VQ(datapath,resultpath):
    """ Calculates and saves conditionwise between-trial variability i.e. average correlation distance between each trial-pairs at each time-point for each conditions separately.

        INPUTS:

        datapath: full path to where the '.mat' file for each subject is stored.
        Each '.mat' file contains array of size: channel x timepoints x trial

        resultpath: full path to store the distances for each subject.

        OUTPUTS:
        VQ : combined fly-by distance for each subject and conditions
        len(VQ) = #subjects,
        Each entry contains a dictionary where key = condition, value = distances
        where, distances =  variance time-series per condition
    """

    # Check if datapath exists
    if(not os.path.isdir(datapath)):
        print('ERROR: Data path does not exist...')
        return

    # Check if result path exists
    if(not os.path.isdir(resultpath)):
        print('ERROR: Result path does not exist \n Create a directory to store results.')
        return

    # List subject files and sort based on names.
    files = os.listdir(datapath)
    files.sort()
    N = len(files)      # Number of Subjects

    VQ = [[] for i in range(N)]  #Between-Trial Variability

    # For each subject, calculate between-trial-vq
    for sub in range(N):
        name = files[sub]
        print('processing subject...' , name)
        data = spio.loadmat(os.path.join(datapath, name))        # Load matrix file

        # List all conditions
        all_vars = spio.whosmat(os.path.join(datapath, name))
        conds = [all_vars[i][0] for i in range(len(all_vars))]

        vq_dist = {}

        # For each experimental condition, calculate VQ

        for args,count in zip(conds,range(len(conds))):
            trials = data[args]                     # List all trials of condition in 'args'
            gfp = np.std(trials,0)                  # Calculate GFP
            trials = np.divide(trials,gfp)      # Divide by GFP
            T = trials.shape[1]

            vq = np.zeros((T))
            for t in range(T):
                vq[t] = np.nanmean(pdist(trials[:,t,:].T,metric='correlation'))
            vq_dist[conds[count]] = vq
        VQ[sub] = vq_dist
        spio.savemat(os.path.join(resultpath, name[:-4] + '_vq.mat'),vq_dist)

    return VQ
