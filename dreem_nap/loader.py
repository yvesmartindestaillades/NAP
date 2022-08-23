from array import array
from typing import Tuple, List, Dict
import pandas as pd
import numpy as np
import pickle
import os
from dreem_nap import manipulator
from tqdm.auto import tqdm

def __load_pickle_to_df(file:str)->pd.DataFrame:
    """Load a pickle file.
    
    Args:
        path (str): the path to the pickle file.
    
    Returns:
        The pickle file content under the dataframe format.    
    """
    assert os.path.exists(file), '{} does not exist.'.format(file)

    with open(file, 'rb') as f:
        mut_hist = pickle.load(f)

    dict_df = {}
    for construct, mh in mut_hist.items():
        dict_df[construct] = __mhs2dict(mh, '_MutationHistogram__bases')

    df = pd.DataFrame.from_dict(dict_df, orient='index').reset_index().drop(columns='index').rename(columns={'name':'construct'})

    return df

def __filter_by_base_cov(df:pd.DataFrame, min_cov_bases:int, samp:str, structure=None, base_type = ['A','C','G','T'], index='all', base_paired=None, flank=None, sub_lib=None):
    """Filter a dataframe by base coverage.
    
    Args:
        df (pd.DataFrame): a dataframe to filter.
        min_cov_bases (int): the minimum base coverage.
    Returns:
        A filtered dataframe.
    """
    args = locals()
    del args['df']
    man = manipulator.Manipulator(df)
    df['worst_cov_bases'] = pd.Series(dtype=float)
    df['min_cov_bases'] = min_cov_bases
    for i, row in tqdm(df.iterrows(), total= len(df.index), unit='construct filtered', colour='green', postfix= 'sample:'+str(samp)):
        if base_type != ['A','C','G','T'] or index != 'all' or base_paired!=None or flank!=None or sub_lib!=None :
            df_loc = man.get_SCC(samp=None, construct=row['construct'], cols=['cov_bases'], cluster=row['cluster'], \
                        **{k:v for k,v in args.items() if k in man.get_SCC.__code__.co_varnames})
            df.loc[i,'worst_cov_bases'] = df_loc['cov_bases'].min()
        else:
            df.loc[i, 'worst_cov_bases'] = np.array(row['cov_bases']).min()
    return df[df['worst_cov_bases'] >= min_cov_bases].reset_index(drop=True)

def __mhs2dict(mhs, drop_attribute:List[str]=[])->dict:
    """Turns the output of DREEM into a 1-level construct-wise index dictionary.

    Args:
        mhs (MutationHistogram): one sample's content under DREEM's MutationHistogram class format. 
        drop_attribute (List[str]): a list of attributes from MutationHistogram class that you don't want into your dictionary
    
    Returns:
        A 1-level dictionary form of the MutationHistogram class.
    """
    mhs_copy = mhs.__dict__.copy()
    for k,v in mhs_copy.items():
        if k in drop_attribute:
            delattr(mhs, k)
        if type(v) == dict:
            for k2,v2 in v.items():
                setattr(mhs, k+'_'+k2, v2)
            delattr(mhs, k)
        if type(v) == np.array:
            setattr(mhs, k, tuple(v))
    return mhs.__dict__

def __filter_by_study(df, samples, name):
    for cons in df.groupby('construct'):
        if len(cons[1]) != len(samples):
            df = df[df['construct'] != cons[0]]
    if df.empty: print('No construct found across all samples for study {}.'.format(name))
    else: print('{} constructs found across all samples for study {}.'.format(len(df.groupby('construct')), name))
    return df

def __add_cols_to_df(df:pd.DataFrame):
    # If no ROI, the entire sequence is considered as ROI.
    if 'cluster' not in df.columns:
        df['cluster'] = 0
    if 'ROI_start' not in df.columns:
        df['ROI_start'] = 0
    if 'ROI_stop' not in df.columns:
        df['ROI_stop'] =  df['sequence'].apply(lambda x: len(x)-1)
    if 'mut_rates' not in df.columns:
        df['mut_rates'] = df.apply(lambda x: np.divide(x['mut_bases'],x['info_bases']), axis=1)
    return df

def __set_indexes_to_0(df:pd.DataFrame):
    for col in df.drop(columns=['sequence','structure']):
        if type(df[col].iloc[0]) in [list, np.ndarray] and len(df[col].iloc[0]) == (len(df['sequence'].iloc[0])+1):
            df[col] = df[col].apply(lambda x: x[1:])
    return df

def df_from_local_files(path_to_data:str, min_cov_bases:int, samples, name, filter_by,  structure, base_type, index, base_paired)->pd.DataFrame:
    args = locals()
    all_df = {}
    assert filter_by in ['study','sample'], 'filter_by must be either study or sample.'
    for s in samples:
        all_df[s] = __load_pickle_to_df(file='{}/{}.p'.format(path_to_data,s))
        all_df[s] = __set_indexes_to_0(all_df[s])
        all_df[s] = __add_cols_to_df(all_df[s])
        all_df[s] = __filter_by_base_cov(all_df[s], min_cov_bases, s, **{k:v for k,v in args.items() if k in manipulator.Manipulator(pd.DataFrame()).get_SCC.__code__.co_varnames})
    df = pd.concat(all_df).reset_index().drop(columns='level_1').rename(columns={'level_0':'samp'})
    if filter_by == 'study':
        df = __filter_by_study(df, samples, name)
    return df

        


