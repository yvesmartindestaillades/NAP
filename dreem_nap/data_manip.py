import pandas as pd
import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
from dreem_nap.study import Study
from dreem_nap import utils




def get_construct_attribute(df:pd.DataFrame, column:str)->pd.DataFrame:    
    """Returns a column of the dataframe w.r.t constructs. 
    
    The value of this column must be only dependent on the construct, and not on the sample. 
    For example, the sequence can be read with this function, but not the mutation rate.  
    
    Args:
        df: a Pandas dataframe.
        attribute: the dataframe's column values that you want to look at.

    Returns:
        A dataframe with the constructs as index and the column as data.
    """   

    return df.set_index('construct').sort_values(column)[column].groupby('construct').apply(lambda x:np.array(x)[0]).sort_values()


def get_roi_info(df:pd.DataFrame, samp:str, construct:int)->pd.DataFrame:
    """Returns a dataframe of the ROI of a specific (samp, construct).

    Args:
        df: a Pandas dataframe.
        samp: a specific sample.
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
    df_use = df.set_index(['samp','construct'])
    start, end = df_use['roi_start_index'].loc[(samp,construct)] , df_use['roi_end_index'].loc[(samp,construct)]     
    mut_per_base = pd.DataFrame({'mut_rate':pd.Series(np.array(df_use[f"mut_bases"].loc[samp, construct][1:])/np.array(df_use[f"info_bases"].loc[samp, construct][1:]), dtype=object),
                            'base':list(df_use['full_sequence'].loc[samp, construct]),
                            'paired': np.array([bool(x != '.') for x in list(df_use['full_structure'].loc[samp,construct])]),\
                            'roi_structure_comparison': pd.Series(list(df_use['roi_structure_comparison'].loc[samp,construct]), index=list(range(start, end)))\
                            ,'roi_deltaG':df_use['roi_deltaG'].loc[samp, construct]})\
                            .dropna()\
                            .reset_index()\
                            .set_index(['base', 'paired', 'roi_structure_comparison','index'])
    return mut_per_base



def columns_to_csv(df:pd.DataFrame, study:Study, columns:List[str], title:str, path:str)->None:
    """Save a subset of a Dataframe to a csv file.

    Args:
        df: a Pandas dataframe.
        study: object containing the samples to save.
        columns: columns to save.
        title: how to name your file.
        path: where to store your file.
    
    Returns:
        The csv file content under the dataframe format.    
    """
    samples = study.samples
    np.seterr(invalid='ignore')
    full_path = utils.make_path(path)
    df_print = df[df.samp.isin(samples)]
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
    all_samples, constructs = list(df.samp.unique()), list(df.construct.unique())
    these_samples, these_constructs = np.array(all_samples)[np.random.randint(0, len(all_samples),n_samples)] , np.array(constructs)[np.random.randint(0, len(constructs), n_constructs)]
    if n_samples == 1:
        these_samples = these_samples[0]
    if n_constructs == 1:
        these_constructs = these_constructs[0]

    return these_samples, these_constructs
