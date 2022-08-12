import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from os.path import exists, dirname
import os, sys
from scipy.stats import linregress
from matplotlib.offsetbox import AnchoredText   

sys.path.append(os.path.abspath(""))
from dreem_nap import data_manip, data_wrangler, util, plot
from typing import Tuple, List



def about_a_sample_construct(df:pd.DataFrame, samp:str, construct:str, path:str)->None:
    """Gives broad information about a given sample-construct, such as the base coverage, mutation histograms, the deltaG plot of the sample and the raw data. 

    Args:
        df (pd.DataFrame): dataframe of interest.
        samp (str): sample of interest.
        construct (int): construct of interest.
        path (str): where you want to store your data.
    """
    # base coverage plot
    plot.base_coverage(df, samp, construct)
    util.save_fig(path=path,
                    title=f"base_cov_{samp}_{construct}")
    plt.close()
    # mut hist plots
    plot.mut_histogram('index',df, samp, construct)
    util.save_fig(path=path,
                    title=f"mut_hist_index_{samp}_{construct}")
    plt.close()
    plot.mut_histogram('partition',df, samp, construct)
    util.save_fig(path=path,
                    title=f"mut_hist_partition_{samp}_{construct}")
    plt.close()
    # deltaG plots
    plot.deltaG(df, samp)
    util.save_fig(path=path,
                    title=f"deltaG_{samp}")
    plt.close()
    # raw data
    data_manip.get_roi_info(df, samp, construct, 
                            bases=['A','C','G','T'],
                            structure='full',
                            roi_range=(-9999, 99999))\
                            .to_csv(f"{path}/info.csv")