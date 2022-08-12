#! usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from os.path import exists, dirname
import os, sys
import numpy as np
import seaborn as sns
import json

path = os.path.dirname('/'.join(os.path.abspath(__file__).split('/')[:-2]))
print(path)
sys.path.append(path)

from dreem_nap import database, util
from dreem_nap.study import Study
import yaml

mpl.use('agg')
mpl.rcParams['figure.dpi'] = 200 # the highest the resolution, the slowest the plotting

####
# SET HYPER PARAMETERS HERE
####

with open('examples/config.yml', 'r') as ymlfile:
    cfg = yaml.safe_load(ymlfile)
for k,v in cfg.items():
    print(k,(30-len(k))*'_',v)

mpl.rcParams['figure.dpi'] = cfg['mpl_rcParams_figure_dpi'] # the highest the resolution, the slowest the plotting
mpl.rcParams["figure.figsize"] = cfg['mpl_rcParams_figure_figsize'] # (width, height) in inches

studies = Study.load_studies(cfg['path_to_studies'])
del studies['all samples']

for study in studies.values():
    study.load_df_from_local_files(path_to_data= cfg['path_to_data'], 
                                   min_cov_bases= cfg['min_cov_bases'])
    for s in study.samples:
        for construct in study.constructs:
            study.mut_histogram(s, construct, 'index')
            util.save_fig(f"/Users/ymdt/src/dreem_nap/data/figs/date/mutation histogram/{study.name}/{s}/{construct}.png")
            plt.close()

