from types import NoneType
import pandas as pd
import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
from dreem_nap import util

class Manipulator():
    def __init__(self, df):
        self._df = df

    ## Building block of Study class
    # data manipulation tools
  
    def assert_structure(self, df, structure):
        assert structure in df.columns, f"Structure {structure} not found"

    def assert_deltaG(self, df, deltaG):
        assert deltaG in df.columns, f"deltaG {deltaG} not found"

    def define_index(self, df, samp, construct, cluster, index):
        if index in ['all','full'] :
            return df.index
        if index == 'roi':
            assert [roi in df.columns for roi in ['ROI_start','ROI_stop']], 'ROI_start and ROI_stop not found'
            return list(range(int(self.get_series(self._df, samp, construct, cluster)['ROI_start']), int(self.get_series(self._df, samp, construct, cluster)['ROI_stop'])))
        if type(index) in [list,tuple]:
            assert [i in list(range(len(self.get_series(self._df, samp, construct, cluster)['sequence']))) for i in index], 'Index out of range'
            return index
        raise ValueError(f"Index {index} not recognized")

    def filter_base_paired(self, df_loc, base_paired):

        if base_paired == True:
            df_loc = df_loc[df_loc['paired'] == True]
        elif base_paired == False:
            df_loc = df_loc[df_loc['paired'] == False]
        elif base_paired == None:
            pass
        return df_loc

    def filter_index(self, df_loc, index):
        return df_loc.loc[index]

    def filter_base_type(self, df_loc, base_type):
        df_loc = pd.concat([df_loc[df_loc['base'] == base] for base in base_type], axis=0)
        return df_loc

    def filter_bases(self, df_loc, base_type, index, base_paired):
        df_loc = self.filter_index(df_loc, index)
        if base_paired != None:
            assert 'paired' in df_loc.columns, 'Give structure to filter on'
            df_loc = self.filter_base_paired(df_loc, base_paired)
        df_loc = self.filter_base_type(df_loc, base_type)
        return df_loc

    def get_SCC(self, samp, construct, cols, cluster=0, structure=None, deltaG=None, base_type = ['A','C','G','T'], index='all', base_paired=None):
        """Returns a dataframe containing the content of a cluster of a sample-construct.

        Args:
            df (pd.Dataframe): A study dataframe.
            samp (str): The sample name.
            construct (str): The construct name.
            cols (list): The columns to be returned.
            cluster (int, optional): The cluster number. Defaults to 0.
            structure (str, optional): Structure to use for the 'paired' column, such as 'structure_ROI_DMS'. Defaults to 'structure'.
            deltaG (str, optional): DeltaG to use for the 'deltaG' column, such as 'deltaG_ens_DMS'. Defaults to 'deltaG'.
            base_type (list, optional): Bases to include. Defaults to ['A','C','G','T'].
            index (str, optional): Index to include. Defaults to 'all'.
            base_paired (_type_, optional): Base pairing to include. None is paired + unpaired, True is paired, False is unpaired. Defaults to None.

        Returns:
            _type_: dataframe containing the content of a cluster of a sample-construct.
        """

        df = self._df.copy()
        
        for cat, input in zip(['structure','deltaG'], [structure,deltaG]):
            if input != None:
                getattr(self, 'assert_'+cat)(df, input)
                assert len([c for c in cols if c.startswith(cat)])==0, f"{input} should be entered as a {cat} argument, not as a col. Cols can't contain {cat}."
                cols += [input]

        for col in cols:
            assert col in df.columns, f"Column {col} not found"

        df_loc = self.get_series(df, samp, construct, cluster)
        for col in [c for c in cols if type(df_loc[c]) in [str]]:
            df_loc[col] = list(df_loc[col])

        df_loc = pd.DataFrame({col: df_loc[col] for col in cols})
        for st in [col for col in cols if col.startswith('structure')]:
            df_loc['paired'] = [{'.':False,'(':True,')':True}[x] for x in df_loc[st]]
            df_loc = df_loc.drop(columns=st)
        
        df_loc = df_loc.rename(columns={'sequence':'base'})
        index = self.define_index(df_loc, samp, construct, cluster, index)
        df_loc = self.filter_bases(df_loc, base_type, index, base_paired)

        return df_loc.sort_index()

    def get_series(self, df, samp, construct, cluster):
        assert len(df_out := df[(df['construct'] == construct)&(df['samp'] == samp)&(df['cluster'] == cluster)]) <= 1, 'More than one row found'
        assert len(df_out) >= 1, 'No row found'
        return df_out.iloc[0]


    def get_construct_attribute(self, column:str)->pd.DataFrame:   

        """Returns a column of the dataframe w.r.t constructs. 
        
        The value of this column must be only dependent on the construct, and not on the sample. 
        For example, the sequence can be read with this function, but not the mutation rate.  
        
        Args:
            attribute: the dataframe's column values that you want to look at.

        Returns:
            A dataframe with the constructs as index and the column as data.
        """   
        if not self._df.empty:
            df = self._df.copy()
            return df.set_index('construct').sort_values(column)[column].groupby('construct').apply(lambda x:np.array(x)[0]).sort_values()

    def filter_df_by_sub_lib(self, sub_lib:str)->pd.DataFrame:
        """Returns a dataframe containing only sublibraries that are contains `sub_lib`.

        Args:
            sub_lib (str): Libraries containing this string will be filtered in.

        Raises:
            Exception: No library contains sub_lib

        Returns:
            pd.DataFrame: A dataframe containing only sublibraries that are contains `sub_lib`.
        """
        df = self._df.copy()

        if sub_lib != None:
            sub_libs = [x for x in df['sub-library'].unique() if sub_lib in x]
            if sub_libs == []:
                raise Exception(f"arg {sub_lib} is not a sub-library of the df")
            df = df[df['sub-library'].isin(sub_libs)]
        return df        

    def columns_to_csv(self, columns:List[str], title:str, path:str)->None:
        """Save a subset of a Dataframe to a csv file.

        Args:
            columns: columns to save.
            title: how to name your file.
            path: where to store your file.
        
        Returns:
            The csv file content under the dataframe format.    
        """
        samples = self.samples
        df = self._df.copy()
        np.seterr(invalid='ignore')
        full_path = util.make_path(path)
        df_print = df[df.samp.isin(samples)]
        df_print = df_print[columns] 
        np.set_printoptions(suppress=True)
        if 'mut_rate' in columns:
            df_print['mut_rate'] = df_print.apply(lambda row: np.float32(np.array(row['mut_bases'])/np.array(row['info_bases'])), axis=1)
        df_print = df_print.drop(columns='index')
        df_print.to_csv(f"{full_path}/{title}")
        return df_print

    def rand_sample_construct(self, n_samples:int=1, n_constructs:int=1)->Tuple[List[str], List[int]]:
        """Pick randomly n_samples samples and n_constructs constructs.

        Args:
            n_samples: the number of samples that you want
            n_construct: the number of constructs that you want. 
        
        Returns:
            Two lists containing the randomly picked elements (resp., samples and constructs).
        """
        df = self._df.copy()
        all_samples = list(df.samp.unique())
        these_samples = np.array(all_samples)[np.random.randint(0, len(all_samples),n_samples)] 
        these_constructs = []
        for samp in these_samples:
            constructs = df[df.samp==samp].construct.unique()
            these_constructs.append(constructs[np.random.randint(0, len(constructs))])
        if n_samples == 1:
            these_samples = these_samples[0]
        if n_constructs == 1:
            these_constructs = these_constructs[0]

        return these_samples, these_constructs    

