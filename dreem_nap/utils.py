import pandas as pd
import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
from dreem_nap.study import Study
import dreem_nap as nap


def roi_range_calc(overlay = 0, roi_range:List[int] = None, roi_bounds:Tuple[int,int] = None, full_bounds:Tuple[int,int] = None):
    """_summary_

    Args:
        overlay (str or int or tuple[int]): extend/shrink roi
            'all': the roi is all bases
            int-type argument: the roi is the subsequence [start_roi_index-overlay, end_roi_index+overlay] 
            tuple[int]-type argument: the roi is the subsequence [start_roi_index+overlay[0], end_roi_index+overlay[1]].
            Defaults to 0.
        roi_range (List[int], optional): Array of base indexes (list[int]). ex: [80, 83, 91]. Defaults to None.
        roi_bounds (Tuple[int,int], optional): Boundaries of the ROI. roi_bounds[0] is included, roi_bounds[1] is excluded. Defaults to None.
        full_bounds (Tuple[int,int], optional): Boundaries of the entire sequence. full_bounds[0] is included, roi_bofull_boundsunds[1] is excluded.. Defaults to None.
    """
    # Select base indexes
    roi_range_name = None
    if overlay == 'all':
        roi_range_name = 'all'
        roi_range = list(range(full_bounds[0], full_bounds[1]))
    if overlay == 0 and roi_range == None:
        roi_range_name = 'roi'
        roi_range = list(range(roi_bounds[0],roi_bounds[1]))
    if overlay != 0 :
        if type(overlay) == int or type(overlay) == float:
            overlay = (-overlay, overlay)
        elif not ((type(overlay) == tuple or type(overlay) == list) and len(overlay) == 2):
            raise f"Unvalid type {type(overlay)} for overlay, please select int."
        roi_range = list(range(roi_bounds[0]+overlay[0], roi_bounds[1]+overlay[1]))
        roi_range_name = f"[start_roi {overlay[0]}, end_roi +{overlay[1]}"
    if roi_range_name == None:
        roi_range_name = 'custom'

    return roi_range, roi_range_name

def make_path(path:str)->str:
    """Create directories until path exists on your computer. Turns the keyword 'date' into today's date.

    Args:
        path: series of directories that you want to create.
    
    Returns:
        Updated path with today's date instead of the keyword 'date'  
    """

    path = os.path.normpath(path)
    path=path.split(os.sep)
    try:
        path[path.index('date')] = str(datetime.datetime.now())[:10]
    except:
        'No date in path'
    full_path = ''
    for repo in path:
        full_path = full_path + f"{repo}/"
        if not exists(full_path):
            os.mkdir(full_path)
    return full_path

def filter_df_by_sub_lib(df:pd.DataFrame, sub_lib:str)->pd.DataFrame:
    """Returns a dataframe containing only sublibraries that are contains `sub_lib`.

    Args:
        df (pd.DataFrame): A standard NAP dataframe.
        sub_lib (str): Libraries containing this string will be filtered in.

    Raises:
        Exception: No library contains sub_lib

    Returns:
        pd.DataFrame: A dataframe containing only sublibraries that are contains `sub_lib`.
    """

    if sub_lib != None:
        sub_libs = [x for x in df['sub-library'].unique() if sub_lib in x]
        if sub_libs == []:
            raise Exception(f"arg {sub_lib} is not a sub-library of the df")
        df = df[df['sub-library'].isin(sub_libs)]
    return df

def gini(x:np.array)->float:
    """Returns Gini index

    Args:
        x (np.array): the array you want the Gini index from

    Returns:
        float: Gini index of the input array
    """
    # (Warning: This is a concise implementation, but it is O(n**2)
    # in time and memory, where n = len(x).  *Don't* pass in huge
    # samples!)

    # Mean absolute difference
    mad = np.abs(np.subtract.outer(x, x)).mean()
    # Relative mean absolute difference
    rmad = mad/np.mean(x)
    # Gini coefficient
    g = 0.5 * rmad
    return g

def save_fig(path:str,title:str, facecolor='white', close=True)->None:
    """Save a matplotlib figure and create the directory if it doesn't exists.

    Args:
        path: where to store your figure.
        title: your figure name.   
        facecolor: color of the background 
    """

    full_path = make_path(path)
    plt.savefig(f"{full_path}/{title}", facecolor=facecolor )
    if close:
        # Clear the current axes.
        plt.cla() 
        # Clear the current figure.
        plt.clf() 
        # Closes all the figure windows.
        plt.close('all')   


def define_figure(title:str, xlabel:str, ylabel:str, figsize:Tuple[float, float])->plt.figure:
    """Define title, labels and size of your figure.

    Args:
        title: matplotlib title
        xlabel: matplotlib xlabel
        ylabel: matplotlib ylabel
        figsize: matplotlib figsize
    """

    fig = plt.figure(figsize=figsize)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return fig
