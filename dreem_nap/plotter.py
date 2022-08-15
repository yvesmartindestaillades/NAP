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
    def __init__(self, df):
        self.__man = manipulator.Manipulator(df)
        self._df = df

    def mut_histogram(self, samp:str, construct:str, cluster:int=0, plot_type:str='index', figsize=(35,7), base_type=['A','C','G','T'], base_index='all', base_paired=None, structure=None, deltaG=None,**kwargs)->None:
        """Plot the mutation rate of a specific (sample, construct).

        Args:
            samp: sample of interest.
            construct: construct of interest.
            plot_type: 'index' or 'partition'. 
                - 'index' uses bases numbers as index and the original construct bases as colors.
                - 'partition' uses original sequence bases as index and the partition of mutated bases as colors.
            figsize: figure size.
            **kwargs: 
                - keyword arguments for matplotlib.pyplot
        
        Returns:
            OutputPlot: output plot data:
                - fig: figure object.
                - ax: axis object.
                - data: plotted data.
        """

        args = locals()
        for attr in ['self','kwargs']:
            del args[attr]
            
        colors = [{'A':'r','C':'b','G':'y','T':'g'}[b] for b in base_type]


        fig = plt.figure(figsize=figsize)

        df_use = self._df.set_index(['samp','construct'])
        
        if not plot_type in ['index','partition']:
            raise Exception(f"{plot_type} must be 'index' or 'partition', please check this argument")

        df_hist = pd.DataFrame()
        df = self.__man.get_SCC(cols = ['sequence','mut_rates']+[f"mod_bases_{base}" for base in base_type],\
                        **{k:v for k,v in args.items() if k in self.__man.get_SCC.__code__.co_varnames})
   
        if plot_type == 'index':  # Plot the mutation rate for each base along the sequence

            mut_per_base = df[['base','mut_rates']].reset_index().set_index(['base', 'index'])
            df_hist.index = df.index

            for base in base_type:
                df_hist[base] = pd.Series(dtype=float)
                df_hist[base] = mut_per_base.loc[base]
            df_hist.dropna(inplace=True, how='all')
            ax = df_hist.plot.bar(stacked=True, color=colors,  figsize=figsize)

        if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
            for base in base_type:
                df_hist[f"mod_bases_{base}"]  = np.array(df_use[f"mod_bases_{base}"].loc[samp, construct][1:])/df_use['info_bases'].loc[samp, construct][1:]

            df_hist.index = list(df_use['sequence'].loc[samp,construct])

            ax = df_hist.plot.bar(stacked=True, color=colors, figsize=figsize)

        plt.title('  '.join([f"{k}: {v}" for k,v in args.items() if \
            k not in ['self','kwargs', 'plot_type','figsize','base_type'] and v is not None]))

        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        
        return OutputPlot(fig, ax, df_hist)