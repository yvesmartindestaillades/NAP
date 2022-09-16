from binascii import a2b_hex
from typing import Tuple, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath(""))
from dreem_nap import manipulator
from itertools import cycle

from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
from tqdm.auto import tqdm


class OutputPlot(object):
    def __init__(self,data, figsize=None, dpi=None) -> None:
        self.fig = plt.figure(figsize=figsize, dpi=dpi)
        self.fig.patch.set_facecolor('white')
        self.ax = plt.axes()
        self.data = data

class Plotter():
    def __init__(self, df):
        self.__man = manipulator.Manipulator(df)

    def mut_histogram(self, samp:str, construct:str, cluster:int=0, plot_type:str='index', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(35,7), **kwargs)->OutputPlot:
        """Plot the mutation rates as histograms.

        Args:
            samp (str): Sample of your sample-construct-cluster.
            construct (str): Construct of your sample-construct-cluster.
            cluster (int, optional): Cluster of your sample-construct-cluster. Defaults to 0. 
            plot_type (str, optional): Colors of the bars.

                * ``'index'`` (default) uses bases numbers as index and the original construct bases as colors.
                * ``'partition'`` uses original sequence bases as index and the partition of mutated bases as colors.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            base_paired (bool, optional): Base-pairing predicition to plot. Defaults to None.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to None.
            show_ci (bool, optional): Show confidence interval on the histogram. Defaults to True.
            figsize (Tuple[int], optional): Figure size. Defaults to (35,7).
            **kwargs: Other arguments to pass to matplotlib.pyplot.

        Raises:
            Exception: plot_type is not ``index`` or ``partition``.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
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
            if show_ci:
                ax = df_hist.plot.bar(y = b, stacked=True, yerr = yerr, legend=b, color=colors[b], ax=ax)
            else:
                ax = df_hist.plot.bar(y = b, stacked=True, legend=b, color=colors[b], ax=ax)

        if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
            df = self.__man.get_SCC(cols = ['sequence','info_bases']+[f"mod_bases_{base}" for base in base_type],\
                        **{k:v for k,v in args.items() if k in self.__man.get_SCC.__code__.co_varnames})

            for base in base_type:
                df_hist[base]  = np.array(df[f"mod_bases_{base}"]/df['info_bases']).astype(float)

            ax = df_hist.plot.bar(stacked=True, color=colors_base, figsize=figsize)

        if len(str(args['index'])) >50:
            args['index'] = str(args['index'])[:50]+'... ]'
        plt.title('  '.join([f"{k}: {v}" for k,v in args.items() if \
            k not in ['self','kwargs', 'plot_type','figsize','base_type'] and not hasattr(plt,k) \
                and v is not None]))

        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        
        return OutputPlot(fig, ax, df)

    def deltaG_sample(self, samp:str, deltaG:str, structure:str, index='all', base_type:List[str]=['A','C','G','T'], flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], figsize:Tuple[int]=(20,5), **kwargs)->OutputPlot:
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            samp (str): Sample of your sample-construct-cluster.
            deltaG (str): DeltaG to use as x axis.
            structure (str, optional): Structure to use for base-paired filtering.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            flank (str, optional): Flank or list of flanks to filter by. Defaults to None.
            sub_lib (str, optional): Sub-library or list of sub-libraries to filter by. Defaults to None.
            max_mutation (float, optional): Maximum mutation rate to plot. Defaults to 0.15.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            figsize (Tuple[int], optional): Figure size. Defaults to (20,5).
            **kwargs: Other arguments to pass to matplotlib.pyplot.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """

        args = locals()
        for attr in ['self','kwargs']:
            del args[attr]
        fit = manipulator.Fit() 
        assert deltaG in self.__man._df.columns, f"deltaG arg {deltaG} isn't in df columns"
        data = self.__man.collect_x_y_paired_unpaired(cols=[deltaG,'mut_rates'], **{k:v for k,v in args.items() if ((k in self.__man.collect_x_y_paired_unpaired.__code__.co_varnames) and (k != 'cols'))})

        out = OutputPlot(data=data, figsize=figsize)       

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

        return out



    def dendrogram(self, samp:str, constructs:str='all', cluster:int=0, metric:str='euclidean', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(10,50), dpi = 600, **kwargs)->OutputPlot:
        """Plot the mutation rates as histograms.

        Args:
            samp (str): Sample(s) of your sample-construct-cluster. A single sample or a list of samples.
            constructs (str): Constructs to plot. Defaults to ``'all'``.
            cluster (int, optional): Cluster of your sample-construct-clusters. Defaults to 0. 
            metric (str, optional): Metric to use for the clustering. Defaults to ``'euclidean'``.
            p (int, optional):
                The ``p`` parameter for ``truncate_mode``.
            truncate_mode (str, optional): 
                The dendrogram can be hard to read when the original
                observation matrix from which the linkage is derived is
                large. Truncation is used to condense the dendrogram. There
                are several modes:

                ``None``
                No truncation is performed (default).
                Note: ``'none'`` is an alias for ``None`` that's kept for
                backward compatibility.

                ``'lastp'``
                The last ``p`` non-singleton clusters formed in the linkage are the
                only non-leaf nodes in the linkage; they correspond to rows
                ``Z[n-p-2:end]`` in ``Z``. All other non-singleton clusters are
                contracted into leaf nodes.

                ``'level'``
                No more than ``p`` levels of the dendrogram tree are displayed.
                A "level" includes all nodes with ``p`` merges from the final merge.

                Note: ``'mtica'`` is an alias for ``'level'`` that's kept for
                backward compatibility.

            color_threshold (double, optional): 
                For brevity, let :math:`t` be the ``color_threshold``.
                Colors all the descendent links below a cluster node
                :math:`k` the same color if :math:`k` is the first node below
                the cut threshold :math:`t`. All links connecting nodes with
                distances greater than or equal to the threshold are colored
                with de default matplotlib color ``'C0'``. If :math:`t` is less
                than or equal to zero, all nodes are colored ``'C0'``.
                If ``color_threshold`` is None or 'default',
                corresponding with MATLAB(TM) behavior, the threshold is set to
                ``0.7*max(Z[:,2])``.

            orientation (str, optional): 
                The direction to plot the dendrogram, which can be any
                of the following strings:

                ``'top'``
                Plots the root at the top, and plot descendent links going downwards.
                (default).

                ``'bottom'``
                Plots the root at the bottom, and plot descendent links going
                upwards.

                ``'left'``
                Plots the root at the left, and plot descendent links going right.

                ``'right'``
                Plots the root at the right, and plot descendent links going left.
            index (optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            base_paired (bool, optional): Base-pairing predicition to plot. Defaults to None.
            structure (str, optio`nal): Structure to use for base_paired filtering. Defaults to None.
            flank (str or list, optional): Flank or list of flanks to filter constructs by. Defaults to None.
            sub_lib (str or list, optional): Sub-library or list of sub-libraries to filter constructs by. Defaults to None.
            figsize (Tuple[int], optional): Figure size. Defaults to (35,7).
            **kwargs: Other arguments to pass to matplotlib.pyplot.

        Returns:
            OutputPlot: The plot object.
        """

        def plot_dendrogram(model, **kwargs):
            """ Create linkage matrix and then plot the dendrogram """

            # create the counts of samples under each node
            counts = np.zeros(model.children_.shape[0])
            n_samples = len(model.labels_)
            for i, merge in enumerate(model.children_):
                current_count = 0
                for child_idx in merge:
                    if child_idx < n_samples:
                        current_count += 1  # leaf node
                    else:
                        current_count += counts[child_idx - n_samples]
                counts[i] = current_count

            linkage_matrix = np.column_stack(
                [model.children_, model.distances_, counts]
            ).astype(float)
            # Plot the corresponding dendrogram
            R = hierarchy.dendrogram(linkage_matrix, labels=model.labels_, **kwargs)
            return R['ivl']


        args = locals()
        del args['self']
        data = self.__man.get_col_across_constructs(col='mut_rates',\
                            **{k:v for k,v in args.items() if k in self.__man.get_SCC.__code__.co_varnames})
        assert metric in ['r2','euclidean','spearman'], f'{metric} is not a valid value for metric'
        affinity = 'euclidean' if metric == 'euclidean' else 'precomputed'
        if affinity == 'precomputed':
            dist_matrix = np.zeros((len(data.index),len(data.index)))
            index = data.index
            data = data.reset_index().drop(columns=['index'])
            for i, u in tqdm(data.iterrows(), total= len(data.index), unit=f' {metric} distance computed', colour='blue'):
                for j, v in data.iterrows():
                    u, v = np.array(u), np.array(v)
                    if metric == 'r2':
                        dist_matrix[i,j] = 1-r2_score(u,v)
                    elif metric == 'spearman':
                        dist_matrix[i,j] = 1-spearmanr(u,v)[0]
            data = pd.DataFrame(dist_matrix, index=index, columns=index)

        out = OutputPlot(data, figsize=figsize)
        # setting distance_threshold=0 ensures we compute the full tree.
        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, affinity=affinity, linkage='ward' if affinity == 'euclidean' else 'complete')
        model = model.fit(data)
        if type(samp) in [str, int]: samp = [str(samp)]
        model.labels_ = [data.index[int(i)] for i in model.labels_]
        plt.title("Hierarchical Clustering Dendrogram")
        # plot the top three levels of the dendrogram
        out.labels = plot_dendrogram(model, truncate_mode="level", orientation='left', **kwargs)
        out.labels.reverse()
        plt.ylabel("Constructs")
        plt.xlabel("Distance")
        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        plt.show()
        return out

