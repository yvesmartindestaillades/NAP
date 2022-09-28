from binascii import a2b_hex
from typing import Tuple, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath(""))
from dreem_nap.manipulator import Manipulator

from dreem_nap.clustering import Clustering
from dreem_nap.util import *
from dreem_nap import deltaG
import plotly

from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from dreem_nap.manipulator import Fit, Manipulator


def mut_histogram(df, samp:str, construct:str, cluster:int=0, structure:str=None, show_ci:bool=True)->OutputPlot:

    mh = Manipulator(df).get_series(df, SubDF.from_locals(locals()))
    xaxis_coordinates = [i for i in range(len(mh.sequence) -1)]

    mut_y = []
    for pos in range(len(mh.sequence)):
        try:
            mut_frac = mh.mut_bases[pos] / mh.info_bases[pos]
        except:
            mut_frac = 0.0
        mut_y.append(mut_frac)
        mut_frac = round(mut_frac, 5)

    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    colors = []
    ref_bases = []
    hover_attr = pd.DataFrame({'mut_rate':mut_y,
                                'base':list(mh.sequence), 
                                'index':list(range(len(mh.sequence))),
                                'paired':[{'.':True, '(':False,')':False}[s] for s in mh.structure]})
    for i in range(len(mh.sequence)):
        if i >= len(mh.sequence)-1:
            continue
        colors.append(cmap[mh.sequence[i - 1]])
        ref_bases.append(mh.sequence[i - 1])
    mut_trace = go.Bar(
            x=xaxis_coordinates,
            y=mut_y,
            text=hover_attr,
            marker=dict(color=colors),
            showlegend=False,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)]),
        )   
    
    if show_ci:
        mut_trace.update(
           error_y=dict(
                type='data',
                symmetric=False,
                array=mh.poisson_high,
                arrayminus=mh.poisson_low
                )
        )

    mut_fig_layout = go.Layout(
            title=f"{mh.samp} - {mh.construct} - {mh.cluster}",
            xaxis=dict(title="Bases"),
            yaxis=dict(title="Mutation rate", range=[0, 0.1]),
            plot_bgcolor="white"

    )
    mut_fig = go.Figure(data=mut_trace, layout=mut_fig_layout)
    seqs = list(mh.sequence)
    if mh.structure is not None:
        db = list(mh.structure)
    else:
        db = " " * len(seqs)
    mut_fig.update_yaxes(
            gridcolor='lightgray',
            linewidth=1,
            linecolor='black',
            mirror=True
    )
    mut_fig.update_xaxes(
            linewidth=1,
            linecolor='black',
            mirror=True
    )
    mut_fig.update_xaxes(
            tickvals=xaxis_coordinates,
            ticktext=["%s<br>%s" % (x, y) for (x, y) in zip(seqs, db)],
            tickangle=0
    )

    plotly.offline.plot(
            mut_fig, filename= "pop_avg.html", auto_open=True,
    )
    return OutputPlot(mh, mut_fig)


class Plotter():
    def __init__(self, df):
        self.__df = df

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
        return 0





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

        return 0

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

        return 0


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

        return 0

