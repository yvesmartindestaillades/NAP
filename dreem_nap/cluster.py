
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

