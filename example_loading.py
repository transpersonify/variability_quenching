from variability import between_trial_VQ, flyby_dist
templatefile = './sample_template/templates_combined.mat'
datapath = './sample_inputs'
resultpath = './sample_outputs'

D,keys = flyby_dist(templatefile,datapath,resultpath)
VQ  = between_trial_VQ(datapath,resultpath)
