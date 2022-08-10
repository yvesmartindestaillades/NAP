from array import array
from typing import Tuple, List, Dict
from dreem.bit_vector import MutationHistogram
import pandas as pd
import numpy as np
import pickle

class Loader:
    def __load_pickle_to_df(self, path:str, samp:str)->pd.DataFrame:
        """Load a pickle file.
        
        Args:
            path (str): the path to the pickle file.
        
        Returns:
            The pickle file content under the dataframe format.    
        """
        with open(path, 'rb') as f:
            mut_hist = pickle.load(f)

        dict_df = {}
        for construct, mh in mut_hist.items():
            dict_df[construct] = self.__mhs2dict(mh, '_MutationHistogram__bases')

        df = pd.DataFrame.from_dict(dict_df, orient='index').reset_index().drop(columns='index').rename(columns={'name':'construct'})
      #  df['samp'] = samp

        return df

    def __filter_by_base_cov(self, df:pd.DataFrame, min_cov_bases:int)->pd.DataFrame:
        """Filter a dataframe by base coverage.
        
        Args:
            df (pd.DataFrame): a dataframe to filter.
            min_cov_bases (int): the minimum base coverage.
        Returns:
            A filtered dataframe.
        """
        df['min_cov_bases'] = df['cov_bases'].apply(lambda x: min(x[1:]))
        return df[df['min_cov_bases'] >= min_cov_bases].reset_index(drop=True)

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

    def __filter_construct(self, df):
        for cons in df.groupby('construct'):
            if len(cons[1]) != len(self.samples):
                df = df[df['construct'] != cons[0]]
        if df.empty: print('No construct found across all samples for study {}.'.format(self.name))
        else: print('{} constructs found across all samples for study {}.'.format(len(df.groupby('construct')), self.name))
        return df

    def load_df_from_local_files(self, path_to_data:str, min_cov_bases:int)->pd.DataFrame:
        all_df = {}
        for s in self.samples:
            all_df[s] = self.__load_pickle_to_df(path='{}/{}/mh.p'.format(path_to_data,s), samp=s)
            all_df[s] = self.__filter_by_base_cov(all_df[s], min_cov_bases)
        self.df = pd.concat(all_df).reset_index().drop(columns='level_1').rename(columns={'level_0':'samp'})
        self.df = self.__filter_construct(self.df)
        return self.df