import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from os.path import exists, dirname
import os, sys
import numpy as np
import seaborn as sns
import json
import yaml
import pickle

path = '/Users/ymdt/src/dreem_nap/'
sys.path.append(path)
from dreem_nap.manipulator import Manipulator
from dreem_nap.study import Study, util

# Config
path_to_data= '/Users/ymdt/src/data/Lauren'
path_to_studies= '/Users/ymdt/src/data/Lauren/studies.csv'
min_cov_bases= 1000
mpl.rcParams['figure.dpi'] = 100 # the highest the resolution, the slowest the plotting
mpl.use('agg')

studies = Study.load_studies(path_to_studies)
study = Study.from_dict(studies['3UTR_v_5UTR'].__dict__)
study.load_df_from_local_files(path_to_data= path_to_data, min_cov_bases = min_cov_bases, filter_by='study')
study.get_df().head()



metric = 'spearman' # can be 'euclidean', 'correlation', 'spearman'
base_type=['A','C']
linkage = 'complete'   
index= list(range(19,42))
        
metric = 'spearman' # can be 'euclidean', 'correlation', 'spearman'
base_type=['A','C']
linkage = 'complete'   
index= list(range(40,50))

#study.plot.cluster_dendrogram(samp=study.samples, index=index, metric=metric, base_type = base_type, linkage=linkage, figsize=(10,70), dpi=300)  
study.plot.deltaG_sample(samp=470,base_type=base_type, deltaG='deltaG_min', structure='structure' )