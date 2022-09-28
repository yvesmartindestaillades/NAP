import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List

from plotly.validators.scatter.marker import SymbolValidator
import plotly.offline as pyo
vals = SymbolValidator().values

def Setshape(x):
        vals = SymbolValidator().values
        return vals[3*x]

class OutputPlot(object):
    def __init__(self,data, fig) -> None:
        self.data = data
        self.fig = fig

class SubDF(object):
    def __init__(self, samp=None,construct=None, cluster=0, structure='structure', base_paired=None, index='all', base_type=['A','C']) -> None:
        for k,v in locals().items():
            setattr(self, k, v)

    def update(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)
            
    @classmethod
    def from_locals(cls, loc):
        return cls(**{k:v for k,v in loc.items() if (k in SubDF.__init__.__code__.co_varnames and k !='self')})


class MplAttr(object):
    def __init__(self, figsize=(10,20), title_fontsize=40, x_ticks_fontsize=30, y_ticks_fontsize=30,  dpi=300) -> None:
        for k,v in locals().items():
            setattr(self, k, v)

    @classmethod
    def from_locals(cls, loc):
        return cls(**{k:v for k,v in loc.items() if (k in MplAttr.__init__.__code__.co_varnames and k !='self')})
        

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

def savefig(file:str, close=True)->None:
    """Save a matplotlib figure and create the directory if it doesn't exists.

    Args:
        file: path+title.
        facecolor: color of the background 
    """

    path = make_path('/'.join(file.split('/')[:-1]))
    plt.savefig(path+file.split('/')[-1])
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



