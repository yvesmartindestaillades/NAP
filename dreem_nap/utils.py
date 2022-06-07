import pandas as pd
import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List

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


def get_construct_attribute(df:pd.DataFrame, attribute:str):   #TODO I don't know wwhat output    
    """Returns columns values that are common to each constructs among the samples - typically, from the RNAstructure output file. 
    
    Args:
        df: a Pandas dataframe.
        column: the dataframe's column values that you want to look at.
    Returns:
        #TODO
    """   

    return df.set_index('construct').sort_values(attribute)[attribute].groupby('construct').apply(lambda x:np.array(x)[0]).sort_values()


def get_roi_info(df:pd.DataFrame, sample:str, construct:int)->pd.DataFrame:
    """Returns a dataframe of the ROI of a specific (sample, construct).

    Args:
        df: a Pandas dataframe.
        sample: a specific sample.
        construct: a specific construct.
    
    Returns:
        Indexes:
            base: A, C, G, T.
            paired: pairing prediction of an RNA structure prediction software.
            roi_structure_comparison: comparison between the pairing prediction of the ROI sequence alone and the entire sequence. 0 means no difference, 1 means difference. 
            index: base 0-index
        Columns:
            mut_rate: mutation rate of this base.
            roi_deltaG: the deltaG of the ROI sequence predicted b a RNA structure prediction software. 
    """

    np.seterr(invalid='ignore')
    df_use = df.set_index(['sample','construct'])
    start, end = df_use['roi_start_index'].loc[(sample,construct)] , df_use['roi_end_index'].loc[(sample,construct)]     
    mut_per_base = pd.DataFrame({'mut_rate':pd.Series(np.array(df_use[f"mut_bases"].loc[sample, construct][1:])/np.array(df_use[f"info_bases"].loc[sample, construct][1:]), dtype=object),
                            'base':list(df_use['full_sequence'].loc[sample, construct]),
                            'paired': np.array([bool(x != '.') for x in list(df_use['full_structure'].loc[sample,construct])]),\
                            'roi_structure_comparison': pd.Series(list(df_use['roi_structure_comparison'].loc[sample,construct]), index=list(range(start, end)))\
                            ,'roi_deltaG':df_use['roi_deltaG'].loc[sample, construct]})\
                            .dropna()\
                            .reset_index()\
                            .set_index(['base', 'paired', 'roi_structure_comparison','index'])
    return mut_per_base

def columns_to_csv(df:pd.DataFrame, samples:List[str], columns:List[str], title:str, path:str)->None:
    """Save a subset of a Dataframe to a csv file.

    Args:
        df: a Pandas dataframe.
        samples: samples to save.
        columns: columns to save.
        title: how to name your file.
        path: where to store your file.
    
    Returns:
        The csv file content under the dataframe format.    
    """

    np.seterr(invalid='ignore')
    full_path = make_path(path)
    df_print = df[df.sample.isin(samples)]
    df_print = df_print[columns] 
    np.set_printoptions(suppress=True)
    if 'mut_rate' in columns:
        df_print['mut_rate'] = df_print.apply(lambda row: np.float32(np.array(row['mut_bases'])/np.array(row['info_bases'])), axis=1)
    df_print = df_print.reset_index().drop(columns='index')
    df_print.to_csv(f"{full_path}/{title}")
    return df_print

def rand_sample_construct(df:pd.DataFrame, n_samples:int=1, n_constructs:int=1)->Tuple[List[str], List[int]]:
    """Pick randomly n_samples samples and n_constructs constructs.

    Args:
        df: a Pandas dataframe to pick samples and constructs from.
        n_samples: the number of samples that you want
        n_construct: the number of constructs that you want. 
    
    Returns:
        Two lists containing the randomly picked elements (resp., samples and constructs).
    """
    all_samples, constructs = list(df.sample.unique()), list(df.construct.unique())
    these_samples, these_constructs = np.array(all_samples)[np.random.randint(0, len(all_samples),n_samples)] , np.array(constructs)[np.random.randint(0, len(constructs), n_constructs)]
    if n_samples == 1:
        these_samples = these_samples[0]
    if n_constructs == 1:
        these_constructs = these_constructs[0]

    return these_samples, these_constructs


def save_fig(path:str,title:str)->None:
    """Save a matplotlib figure and create the directory if it doesn't exists.

    Args:
        path: where to store your figure.
        title: your figure name.    
    """

    full_path = make_path(path)
    plt.savefig(f"{full_path}/{title}")


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