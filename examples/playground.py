
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, sys
import numpy as np
from tqdm.auto import tqdm

sys.path.append('/Users/ymdt/src/dreem_nap/')
from dreem_nap import manipulator 
from dreem_nap.study import Study, util
from dreem_nap.util import *
from itertools import cycle
import plotly.express as px
from scipy.optimize import curve_fit
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from dreem_nap.manipulator import Fit

# Create a study
samples_csv = pd.read_csv('~/src/data/Jordan/samples.csv')

dms = Study.from_dict({'name': 'dms',
                         'description': 'Change the DMS concentration', 
                         'samples': list(samples_csv['sample']), 
                         'label': 'Na quantity [M]', 
                         'conditions': list(samples_csv['DMS_conc_mM'])})

# Load data
dms.load_df_from_local_files(path_to_data='/Users/ymdt/src/data/Jordan', 
                              min_cov_bases= 1000, 
                              filter_by='sample')


    
dms._df.head(3)


def per_base(df:pd.DataFrame, construct:str, experimental_variable:str, structure:str='structure', index='all', base_type=['A','C','G','T'], flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], figsize:Tuple[int]=(20,5), title_fontsize=40, xticks_fontsize=30, yticks_fontsize=30, **kwargs)->OutputPlot:

    fit = manipulator.Fit()
    man = manipulator.Manipulator(df)
    data = pd.DataFrame()
    sub_df= SubDF.from_locals(locals())
    assert experimental_variable in df.columns, f"{experimental_variable=} isn't in the study columns"
    hover_attr = ['sequence','mut_rates','samp',experimental_variable,]
    for row in df[df.construct==construct].itertuples():
        sub_df.update(construct = row.construct, cluster=row.cluster)
        data = pd.concat((data,man.get_SCC(cols=hover_attr.copy(),sub_df=sub_df,can_be_empty=True).reset_index()))
    return data


print(per_base(dms._df, '410-O-flank_1=hp14', 'DMS_conc_mM'))