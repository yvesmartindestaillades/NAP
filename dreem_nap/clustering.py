
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from typing import Tuple, List
from dreem_nap.plotter import OutputPlot
from dreem_nap.manipulator import Manipulator

def __plot_dendrogram(model, **kwargs):
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
    hierarchy.dendrogram(linkage_matrix, labels=model.labels_, **kwargs)
    


class Clustering:
    def __init__(self, df) -> None:
        self.__man = Manipulator(df)

    def dendrogram(self, samp:str, constructs:str='all', cluster:int=0, metric:str='euclidean', index='all', base_type:List[str]=['A','C','G','T'], base_paired:bool=None, structure:str=None, show_ci:bool=True, figsize:Tuple[int]=(10,50), dpi = 600, **kwargs)->OutputPlot:
        args = locals()
        del args['self']
        data = self.__man.get_col_across_constructs(col='mut_rates',\
                            **{k:v for k,v in args.items() if k in self.__man.get_SCC.__code__.co_varnames})

        # setting distance_threshold=0 ensures we compute the full tree.
        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)
        model = model.fit(data)
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.axes()
        if type(samp) in [str, int]: samp = [str(samp)]
        model.labels_ = [data.index[int(i)] for i in model.labels_]
        plt.title("Hierarchical Clustering Dendrogram")
        # plot the top three levels of the dendrogram
        __plot_dendrogram(model, truncate_mode="level", orientation='left', **kwargs)

        plt.ylabel("Constructs")
        plt.xlabel("Distance")
        plt.show()

        out = OutputPlot(fig, ax, data)
        out.labels = model.labels_
        return out