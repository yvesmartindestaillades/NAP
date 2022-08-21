import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath(""))
from dreem_nap import manipulator
from itertools import cycle

class OutputPlot(object):
    def __init__(self, fig, ax, data) -> None:
        self.fig = fig
        self.ax = ax
        self.data = data

class Plotter():
    def __init__(self, df):
        self.__man = manipulator.Manipulator(df)

    def mut_histogram(self, samp:str, construct:str, cluster:int=0, plot_type:str='index', figsize=(35,7), base_type=['A','C','G','T'], index='all', base_paired=None, structure=None, **kwargs)->None:
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
        
        if 'facecolor' not in kwargs:
            kwargs['facecolor'] = 'w'

        colors = {'A':'r','C':'b','G':'y','T':'g'}
        colors_base = [colors[b] for b in base_type]

        fig = plt.figure(figsize=figsize)
        ax = plt.axes()

        if not plot_type in ['index','partition']:
            raise Exception(f"{plot_type} must be 'index' or 'partition', please check this argument")

        df_hist = pd.DataFrame()
        
        df = self.__man.get_SCC(cols = ['sequence','mut_rates','poisson_low','poisson_high'],\
                        **{k:v for k,v in args.items() if k in self.__man.get_SCC.__code__.co_varnames})

        mut_per_base = df[['base','mut_rates']].reset_index().set_index(['base', 'index'])
        df_err_low = df[['base','poisson_low']].reset_index().set_index(['base', 'index'])
        df_err_high = df[['base','poisson_high']].reset_index().set_index(['base', 'index'])
        df_hist.index = df.index

        for base in base_type:
            df_hist[base], df_hist[base+'_min'], df_hist[base+'_max'] = pd.Series(dtype=float), pd.Series(dtype=float), pd.Series(dtype=float)
            try:
                df_hist[base] = mut_per_base.loc[base]
                df_hist[base+'_min'] =  df_err_low.loc[base] #+ mut_per_base.loc[base] 
                df_hist[base+'_max'] =  df_err_high.loc[base] #- mut_per_base.loc[base]
            except:
                f"No mutations for base {base}"
        df_hist.dropna(inplace=True, how='all')
        for b in base_type:
            yerr = np.array(df_hist[[b+'_min',b+'_max']]).T
            ax = df_hist.plot.bar(y = b, stacked=True, yerr = yerr, legend=b, color=colors[b], ax=ax)

        if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
            df = self.__man.get_SCC(cols = ['sequence','info_bases']+[f"mod_bases_{base}" for base in base_type],\
                        **{k:v for k,v in args.items() if k in self.__man.get_SCC.__code__.co_varnames})

            for base in base_type:
                df_hist[base]  = np.array(df[f"mod_bases_{base}"]/df['info_bases']).astype(float)

            ax = df_hist.plot.bar(stacked=True,
             color=colors_base, figsize=figsize)

        if len(str(args['index'])) >50:
            args['index'] = str(args['index'])[:50]+'... ]'
        plt.title('  '.join([f"{k}: {v}" for k,v in args.items() if \
            k not in ['self','kwargs', 'plot_type','figsize','base_type'] and not hasattr(plt,k) \
                and v is not None]))

        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        
        return OutputPlot(fig, ax, df_hist)


    def deltaG_sample(self, samp:str, structure, deltaG, figsize=(20,5), base_type=['A','C','G','T'], index='all', flank=None, sub_lib=None, max_mutation= 0.15,models=[],**kwargs):
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.
        
        Args:
        """
        args = locals()
        for attr in ['self','kwargs']:
            del args[attr]
        fit = manipulator.Fit()        
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
        assert deltaG in self.__man._df.columns, f"deltaG arg {deltaG} isn't in df columns"
        data = self.__man.collect_x_y_paired_unpaired(cols=[deltaG,'mut_rates'], **{k:v for k,v in args.items() if ((k in self.__man.collect_x_y_paired_unpaired.__code__.co_varnames) and (k != 'cols'))})

        for is_paired, color in zip([True,False],['g','r']):
            plt.plot(data[is_paired]['x'],data[is_paired]['y'], color+'.')

        for color, is_paired, prefix in zip(['g','r'],[True,False], ['Paired bases ','Unpaired bases ']):
            style = cycle(['-','--',':','-.'])
            for m in models:
                plt.plot(*fit.predict(data[is_paired]['x'],data[is_paired]['y'],  m, prefix=prefix), color=color, linestyle=next(style))
            
        plt.legend(['Paired bases data','Unpaired bases data'] + fit.legend )

        if len(str(args['index'])) >50:
            args['index'] = str(args['index'])[:50]+'... ]'
        plt.title('  '.join([f"{k}: {v}" for k,v in args.items() if \
            k not in ['self','kwargs', 'plot_type','figsize','base_type','df']\
            and not hasattr(plt,k) and v is not None]))
        plt.xlabel(f'DeltaG: {deltaG} [Kcal/mol]')
        plt.ylabel(f'Mutation rate')

        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 

        return OutputPlot(fig, ax, data)


