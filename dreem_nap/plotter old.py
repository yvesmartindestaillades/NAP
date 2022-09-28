from binascii import a2b_hex
from typing import Tuple, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath(""))
from dreem_nap import manipulator

from dreem_nap.clustering import Clustering
from dreem_nap.util import *
from dreem_nap.deltaG import DeltaG

class Plotter():
    def __init__(self, df):
        self.__man = manipulator.Manipulator(df)

    def mut_histogram(self, samp:str, construct:str, cluster:int=0, plot_type:str='index', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(35,7), title_fontsize=40, xticks_fontsize=10, yticks_fontsize=30, **kwargs)->OutputPlot:
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
            title_fontsize (int, optional): Title font size. Defaults to 40.
            yticks_fontsize (int, optional): Ytick font size. Defaults to 30.
            **kwargs: Other arguments to pass to matplotlib.pyplot.

        Raises:
            Exception: plot_type is not ``index`` or ``partition``.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """
        sub_df = SubDF.from_locals(locals())
        mpl_attr = MplAttr.from_locals(locals())
        
        if 'facecolor' not in kwargs:
            kwargs['facecolor'] = 'w'

        colors = {'A':'r','C':'b','G':'y','T':'g'}
        colors_base = [colors[b] for b in base_type]

        if not plot_type in ['index','partition']:
            raise Exception(f"{plot_type} must be 'index' or 'partition', please check this argument")

        df_hist = pd.DataFrame()
        
        df = self.__man.get_SCC(cols = ['sequence','mut_rates','poisson_low','poisson_high'],sub_df=sub_df)
        out = OutputPlot(df, mpl_attr)
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
                out.ax = df_hist.plot.bar(y = b, stacked=True, yerr = yerr, legend=b, color=colors[b], ax=out.ax)
            else:
                out.ax = df_hist.plot.bar(y = b, stacked=True, legend=b, color=colors[b], ax=out.ax)

        if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
            df = self.__man.get_SCC(cols = ['sequence','info_bases']+[f"mod_bases_{base}" for base in base_type],\
                        sub_df=sub_df)

            for base in base_type:
                df_hist[base]  = np.array(df[f"mod_bases_{base}"]/df['info_bases']).astype(float)

            out.ax = df_hist.plot.bar(stacked=True, color=colors_base, figsize=figsize, ax=out.ax)
        

        if len(str(index)) >50:
            index = str(index)[:50]+'... ]'
        plt.title('  '.join([f"{k}: {v}" for k,v in sub_df.__dict__.items() if \
            k not in ['base_type','self'] and not hasattr(plt,k) \
                and v is not None]),
                fontsize=title_fontsize)
        plt.xticks(fontsize=xticks_fontsize)
        plt.yticks(fontsize=yticks_fontsize)
        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        
        return out

    def deltaG_sample(self, samp:str, deltaG:str, structure:str, index='all', base_type:List[str]=['A','C','G','T'], flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], figsize:Tuple[int]=(20,5), title_fontsize=40, xticks_fontsize=30, yticks_fontsize=30, **kwargs)->OutputPlot:
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
            title_fontsize (int, optional): Title font size. Defaults to 40.
            xyticks_fontsize (int, optional): X and Y ticks font size. Defaults to 30.
            **kwargs: Other arguments to pass to matplotlib.pyplot.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """
        sub_df = SubDF.from_locals(locals())
        mpl_attr = MplAttr.from_locals(locals())

        return DeltaG(self.__man).per_sample(sub_df, mpl_attr, deltaG, flank, sub_lib, max_mutation, models, **kwargs)

    def deltaG_basewise(self, samp:str, construct:str, cluster:str, deltaG:str, structure:str, index='all', base_type:List[str]=['A','C','G','T'], max_mutation:float= 0.15, models:List[str]=[], figsize:Tuple[int]=(20,5), title_fontsize=40, xticks_fontsize=30, yticks_fontsize=30, **kwargs)->OutputPlot:
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            samp (str): Sample of your sample-construct-cluster.
            construct (str): Construct of your sample-construct-cluster.
            cluster (str): Cluster of your sample-construct-cluster.
            deltaG (str): DeltaG to use as x axis.
            structure (str, optional): Structure to use for base-paired filtering.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            flank (str, optional): Flank or list of flanks to filter by. Defaults to None.
            sub_lib (str, optional): Sub-library or list of sub-libraries to filter by. Defaults to None.
            max_mutation (float, optional): Maximum mutation rate to plot. Defaults to 0.15.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            figsize (Tuple[int], optional): Figure size. Defaults to (20,5).
            title_fontsize (int, optional): Title font size. Defaults to 40.
            xyticks_fontsize (int, optional): X and Y ticks font size. Defaults to 30.
            **kwargs: Other arguments to pass to matplotlib.pyplot.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """
        sub_df = SubDF.from_locals(locals())
        mpl_attr = MplAttr.from_locals(locals())

        return DeltaG(self.__man).per_base(sub_df, mpl_attr, deltaG, max_mutation, models, **kwargs)



    def cluster_dendrogram(self, samp:str, construct:str='all', cluster:int=0, metric:str='euclidean', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(10,50), dpi = 600, distance_threshold=0, n_clusters=None, **kwargs)->OutputPlot:
        """Plot the mutation rates as histograms.

        Args:
            samp (str): Sample(s) of your sample-construct-cluster. A single sample or a list of samples.
            constructs (str): Constructs to plot. Defaults to ``'all'``.
            cluster (int, optional): Cluster of your sample-construct-clusters. Defaults to 0. 
            metric (str, optional): can be  'braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
        'jaccard', 'jensenshannon', 'kulsinski', 'kulczynski1',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
        'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
        'spearman', 'sqeuclidean', 'yule' . Defaults to 'euclidean'.
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

        sub_df = SubDF.from_locals(locals())
        mpl_attr = MplAttr.from_locals(locals())
    
        return Clustering(self.__man).dendrogram(sub_df, mpl_attr, metric=metric, distance_threshold=distance_threshold, n_clusters=n_clusters, **kwargs)

    def cluster_distance_matrix_heatmap(self, samp:str, construct:str='all', cluster:int=0, metric:str='euclidean', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(10,10), dpi = 600, distance_threshold=0, n_clusters=None, **kwargs)->OutputPlot:
        """Plot the distance matrix heatmap.

        Args:
            samp (str): samples to use
            construct (str, optional): constructs to use. Defaults to 'all'.
            cluster (int, optional): clusters to use. Defaults to 0.
            metric (str, optional): can be  'braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
        'jaccard', 'jensenshannon', 'kulsinski', 'kulczynski1',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
        'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
        'spearman', 'sqeuclidean', 'yule' . Defaults to 'euclidean'.
            index (str, optional): Defaults to 'all'.
            base_type (List[str], optional): Defaults to ['A','C','G','T'].
            base_paired (bool, optional): Defaults to None.
            structure (str, optional): Defaults to None.
            figsize (Tuple[int], optional): Defaults to (15,15).
            dpi (int, optional): Defaults to 300.
            distance_threshold (int, optional): float, default=None
        The linkage distance threshold above which, clusters will not be
        merged. If not ``None``, ``n_clusters`` must be ``None`` and
        ``compute_full_tree`` must be ``True``. Defaults to 0.
            n_clusters (int, optional): Number of clusters to find. Defaults to None.

        Returns:
            OutputPlot: _description_
        """

        sub_df = SubDF.from_locals(locals())
        mpl_attr = MplAttr.from_locals(locals())
        
        return Clustering(self.__man).distance_matrix_heatmap( sub_df, mpl_attr, metric=metric, distance_threshold=distance_threshold, n_clusters=n_clusters,**kwargs)

    def cluster_mut_rates_heatmap(self, samp:str, construct:str='all', cluster:int=0, metric:str='euclidean', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(10,10), dpi = 600, distance_threshold=0, n_clusters=None, **kwargs)->OutputPlot:
        """Plot the mutation rate heatmap ordered by a hierarchical cluster.

        Args:
            samp (str): samples to use
            construct (str, optional): constructs to use. Defaults to 'all'.
            cluster (int, optional): clusters to use. Defaults to 0.
            metric (str, optional): can be  'braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
        'jaccard', 'jensenshannon', 'kulsinski', 'kulczynski1',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
        'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
        'spearman', 'sqeuclidean', 'yule' . Defaults to 'euclidean'.
            index (str, optional): Defaults to 'all'.
            base_type (List[str], optional): Defaults to ['A','C','G','T'].
            base_paired (bool, optional): Defaults to None.
            structure (str, optional): Defaults to None.
            figsize (Tuple[int], optional): Defaults to (15,15).
            dpi (int, optional): Defaults to 300.
            distance_threshold (int, optional): float, default=None
        The linkage distance threshold above which, clusters will not be
        merged. If not ``None``, ``n_clusters`` must be ``None`` and
        ``compute_full_tree`` must be ``True``. Defaults to 0.
            n_clusters (int, optional): Number of clusters to find. Defaults to None.

        Returns:
            OutputPlot: _description_
        """
        sub_df = SubDF.from_locals(locals())
        mpl_attr = MplAttr.from_locals(locals())

        return Clustering(self.__man).mut_rates_heatmap( sub_df, mpl_attr,metric=metric, distance_threshold=distance_threshold, n_clusters=n_clusters, **kwargs)

