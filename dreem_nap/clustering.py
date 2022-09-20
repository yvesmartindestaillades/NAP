
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from typing import Tuple, List
from dreem_nap.util import *
from dreem_nap.manipulator import Manipulator
from scipy.stats import spearmanr
from tqdm.auto import tqdm
import scipy
import seaborn as sns

class Clustering:
    def __init__(self, man) -> None:
        self.__man = man

    def mut_rates_heatmap(self,  sub_df:SubDF, mpl_attr:MplAttr, metric:str='euclidean', distance_threshold=0, n_clusters=None, linkage='complete', **kwargs)->OutputPlot:
        args = locals()
        del args['self']
    
        data = self.__man.get_col_across_constructs(col='mut_rates',sub_df=sub_df)
        model = self.__build_model(data, metric, distance_threshold, n_clusters, linkage, **kwargs)

        labels = self.__get_labels(model)
        out = OutputPlot(data,mpl_attr=mpl_attr)
        out.labels = labels
        out.data = data.reindex(out.labels)
        plt.title("Mutation rates sorted with hierarchical clustering by distance ".format(metric))
        plt.xlabel("Bases")
        plt.ylabel("Constructs")
        sns.heatmap(out.data, cmap='viridis', **{k:v for k,v in kwargs.items() if k in sns.heatmap.__code__.co_varnames})
        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        plt.tight_layout()
        return out

    def distance_matrix_heatmap(self, sub_df:SubDF, mpl_attr:MplAttr, metric:str='euclidean', distance_threshold=0, n_clusters=None, linkage='complete', **kwargs)->OutputPlot:
        args = locals()
        del args['self']    
        data = self.__man.get_col_across_constructs(col='mut_rates',sub_df=sub_df)
        dist_matrix = self.__compute_matrix(data, metric=metric)
        model = self.__build_model(data, metric, distance_threshold, n_clusters, linkage, **kwargs)
        labels = self.__get_labels(model)

        out = OutputPlot(data,mpl_attr=mpl_attr)
        out.labels = labels
        for axis in [0,1]:
            dist_matrix = dist_matrix.reindex(out.labels, axis=axis)
        out.data = data.reindex(out.labels)
        plt.title("Distance Matrix {}".format(metric))
        plt.xlabel("Constructs")
        plt.ylabel("Constructs")
        sns.heatmap(dist_matrix, cmap='viridis', **{k:v for k,v in kwargs.items() if k in sns.heatmap.__code__.co_varnames})
        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        plt.tight_layout()
        return out

    def dendrogram(self, sub_df:SubDF, mpl_attr:MplAttr, metric:str='euclidean', distance_threshold=0, n_clusters=None, linkage='complete', **kwargs)->OutputPlot:
        args = locals()
        del args['self']
        data = self.__man.get_col_across_constructs(col='mut_rates', sub_df=sub_df)
        model = self.__build_model(data, metric, distance_threshold, n_clusters, linkage, **kwargs)
        out = self.__plot_dendrogram(model, data, mpl_attr=mpl_attr,**kwargs)
        out.data = out.data.reindex(out.labels)
        plt.title("Hierarchical Clustering Dendrogram")
        plt.ylabel("Constructs")
        plt.xlabel("Distance {}".format(metric))
        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        plt.tight_layout()
        return out

    def __plot_dendrogram(self, model, data=None, mpl_attr=None, close_plot=False, **kwargs):
        """ Create linkage matrix and then plot the dendrogram """
        # create the counts of samples under each node
        out = OutputPlot(data=data,mpl_attr=mpl_attr)
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
        dend = hierarchy.dendrogram(linkage_matrix, labels=model.labels_, truncate_mode="level", orientation='left',**{k:v for k,v in kwargs.items() if k in hierarchy.dendrogram.__code__.co_varnames})
        if close_plot: 
            plt.close()
        dend['ivl'].reverse()
        out.labels = dend['ivl']
        return out

    def __compute_matrix(self, data, metric):
        """Compute the distance matrix between the data.

        Args:
            data (np.array): the data you want to compute the distance matrix from
            metric (str): the metric you want to use to compute the distance matrix

        Returns:
            np.array: the distance matrix
        """
        if metric != 'spearman':
            dist_matrix = scipy.spatial.distance.squareform(hierarchy.distance.pdist(data, metric=metric))
        else:
            dist_matrix = np.zeros((len(data.index),len(data.index)))
            data_iter = data.reset_index().drop(columns='index')
            for i, u in tqdm(data_iter.iterrows(), total= len(data.index), unit=f' {metric} distance computed', colour='blue'):
                for j, v in data_iter.iterrows():
                    dist_matrix[i,j] = 1 - spearmanr(np.array(u), np.array(v))[0]

        return pd.DataFrame(dist_matrix, index=data.index, columns=data.index)

    def __build_model(self, data, metric, distance_threshold, n_clusters, linkage, **kwargs):
        """Build a model of the data.

        Args:
            data (np.array): the data you want to cluster
            metric (str): the metric you want to use to compute the distance matrix
            **kwargs: keyword arguments for the AgglomerativeClustering algorithm

        Returns:
            model: the model of the data
        """
        model = AgglomerativeClustering(distance_threshold=distance_threshold, 
                                        n_clusters=n_clusters, 
                                        affinity='precomputed',
                                        linkage=linkage,
                                        **{k:v for k,v in kwargs.items() if k in AgglomerativeClustering.__init__.__code__.co_varnames})
        model = model.fit(self.__compute_matrix(data, metric=metric))
        model.labels_ = [model.feature_names_in_[i] for i in model.labels_]
        return model

    def __get_labels(self, model):
        """Sort the labels of the clusters.

        Args:
            labels (np.array): the labels of the clusters
            data (np.array): the data you want to sort by

        Returns:
            np.array: the sorted labels
        """
        out = self.__plot_dendrogram(model, close_plot=True, figsize=(2,2))
        return out.labels