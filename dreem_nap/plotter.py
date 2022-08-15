from turtle import filling
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
from dreem_nap import manipulator, util
from typing import Tuple, List
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

class OutputPlot(object):
    def __init__(self, fig, ax, data) -> None:
        self.fig = fig
        self.ax = ax
        self.data = data

class Plotter():
    def mut_histogram(self, samp:str, construct:str, plot_type:str='index', figsize=(35,7), **kwargs)->None:
        """Plot the mutation rate of a specific (sample, construct).

        Args:
            samp: sample of interest.
            construct: construct of interest.
            plot_type: 'index' or 'partition'. 
                - 'index' uses bases numbers as index and the original construct bases as colors.
                - 'partition' uses original sequence bases as index and the partition of mutated bases as colors.
            figsize: figure size.
            **kwargs: keyword arguments for matplotlib.pyplot
        
        Returns:
            OutputPlot: output plot data:
                - fig: figure object.
                - ax: axis object.
                - data: plotted data.
        """
        fig = plt.figure(figsize=figsize)

        df_use = self.df.set_index(['samp','construct'])
        
        if not plot_type in ['index','partition']:
            raise Exception(f"{plot_type} must be 'index' or 'partition', please check this argument")

        df_hist = pd.DataFrame()

        if plot_type == 'index':  # Plot the mutation rate for each base along the sequence

            mut_per_base = pd.DataFrame({'mut_rates': df_use['mut_rates'].loc[samp, construct]
                                        ,'base':list(df_use['sequence'].loc[samp, construct])})\
                                        .reset_index()\
                                        .set_index(['base', 'index'])
            df_hist.index = mut_per_base.reset_index()['index']

            for base in ['A','C','G','T']:
                df_hist[base] = pd.Series(dtype=float)
                df_hist[base] = mut_per_base.loc[base]

            ax = df_hist.plot.bar(stacked=True, color=['r','b','y','g'],  figsize=figsize)
            plt.title(f"sample {samp}, construct {construct}")

        if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
            for base in ['A','C','G','T']:
                df_hist[f"mod_bases_{base}"]  = np.array(df_use[f"mod_bases_{base}"].loc[samp, construct][1:])/df_use['info_bases'].loc[samp, construct][1:]

            df_hist.index = list(df_use['sequence'].loc[samp,construct])

            ax = df_hist.plot.bar(stacked=True, color=['r','b','y','g'], figsize=figsize)
        
        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 

        return OutputPlot(fig, ax, df_hist)