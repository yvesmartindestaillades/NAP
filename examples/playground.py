#! usr/bin/env python3

import struct
from zlib import DEF_BUF_SIZE
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from os.path import exists, dirname
import os, sys
import numpy as np
import seaborn as sns
import json

path = '/Users/ymdt/src/dreem_nap/'

sys.path.append(path)

from dreem_nap.study import Study, util
import yaml
import pickle


with open(path+'config.yml', 'r') as ymlfile:
    cfg = yaml.safe_load(ymlfile)

mpl.rcParams['figure.dpi'] = cfg['mpl_rcParams_figure_dpi'] # the highest the resolution, the slowest the plotting

####
# SET HYPER PARAMETERS HERE
####

studies = Study.load_studies(cfg['path_to_studies'])
study = Study().from_dict(studies['temperature'].__dict__)

with open(path+'data/temperature_df.p','rb') as f:
    study.set_df(pickle.load(f))
    f.close()


samp, construct = 'B8','12419'


study._df['cluster'] = 0
study.constructs = study._df['construct'].unique()

# base_type = ['A','C','G','T']
# base_index = 'roi', 'all', [93,95,96]
# base_paired = True, False or None (=both) # default is None 
# figsize = (25, 7) # custom by plot type

# structure = "structure_ROI"

# deltaG = "deltaG_ens_DMS"

# cluster = 0, 1, 2


# possible indexes
# roi, all, [93,95,96]

#print(dir(study.plot))
study.plot.mut_histogram(samp=samp, construct=construct, \
    plot_type='index',\
    base_type=['A','G','C','T'], index='all', base_paired=None,\
    structure = 'structure', figsize=(25, 7))
