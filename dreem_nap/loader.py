from array import array
from typing import Tuple, List, Dict
from dreem.bit_vector import MutationHistogram
import pandas as pd
import numpy as np
import pickle
import os

class Loader:
    def __load_pickle_to_df(self, file:str)->pd.DataFrame:
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
            dict_df[construct] = self.__mhs2dict(mh, '_MutationHistogram__bases')

        df = pd.DataFrame.from_dict(dict_df, orient='index').reset_index().drop(columns='index').rename(columns={'name':'construct'})

        return df

    def __filter_by_base_cov(self, df:pd.DataFrame, min_cov_bases:int)->pd.DataFrame:
        """Filter a dataframe by base coverage.
        
        Args:
            df (pd.DataFrame): a dataframe to filter.
            min_cov_bases (int): the minimum base coverage.
        Returns:
            A filtered dataframe.
        """
        df['min_cov_bases'] = min_cov_bases
        df['cov_bases_roi'] = df.apply(lambda x: min(x['cov_bases'][int(x['ROI_start']):int(x['ROI_stop'])]), axis=1)
        return df[df['cov_bases_roi'] >= min_cov_bases].reset_index(drop=True)

    def __mhs2dict(self, mhs:MutationHistogram, drop_attribute:List[str]=[])->dict:
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

    def __filter_by_study(self, df, samples, name):
        for cons in df.groupby('construct'):
            if len(cons[1]) != len(samples):
                df = df[df['construct'] != cons[0]]
        if df.empty: print('No construct found across all samples for study {}.'.format(name))
        else: print('{} constructs found across all samples for study {}.'.format(len(df.groupby('construct')), name))
        return df

    def __add_cols_to_df(self, df:pd.DataFrame):
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

    def __set_indexes_to_0(self, df:pd.DataFrame):
        for col in df.drop(columns=['sequence','structure']):
            if type(df[col].iloc[0]) in [list, np.ndarray] and len(df[col].iloc[0]) == (len(df['sequence'].iloc[0])+1):
                df[col] = df[col].apply(lambda x: x[1:])
        return df

    def df_from_local_files(self, path_to_data:str, min_cov_bases:int, samples, name, filter_by)->pd.DataFrame:
        all_df = {}
        assert filter_by in ['study','sample'], 'filter_by must be either study or sample.'
        for s in samples:
            all_df[s] = self.__load_pickle_to_df(file='{}/{}.p'.format(path_to_data,s))
            all_df[s] = self.__set_indexes_to_0(all_df[s])
            all_df[s] = self.__add_cols_to_df(all_df[s])
            all_df[s] = self.__filter_by_base_cov(all_df[s], min_cov_bases)
        df = pd.concat(all_df).reset_index().drop(columns='level_1').rename(columns={'level_0':'samp'})
        if filter_by == 'study':
            df = self.__filter_by_study(df, samples, name)
        return df

        


