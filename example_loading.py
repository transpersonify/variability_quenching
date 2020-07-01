from variability import flyby_dist 
templatefile = '~/Desktop/Backup_05_11/Milad_data/templates/templates_combined.mat'
datapath = '/home/sn254804/Desktop/Backup_05_11/Milad_data/datafiles/'
resultpath = '/home/sn254804/Desktop/Backup_05_11/Milad_data/flybys/'

D,keys = flyby_dist(templatefile,datapath,resultpath)
