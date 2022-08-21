import pandas as pd
import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
from dreem_nap import util

MAX_MUTATION = 0.3

from typing import Tuple
from scipy.optimize import curve_fit
from itertools import cycle

class Fit(object):
    def __init__(self) -> None:
        self.legend = []
    
    def clean_legend(self):
        self.legend = []

    def predict(self, x, y, model, prefix='', suffix=''):
        fit = self.fit(x,y,model, prefix, suffix)
        m = eval(model)
        return np.sort(x), m(np.sort(x),*fit)

    def fit(self, x,y, model, prefix, suffix):
        fit = curve_fit(eval(model), x, y)[0]
        self._generate_legend(fit, model, prefix, suffix)
        return fit

    def _generate_legend(self, fit, m, prefix, suffix):
        slice_m = lambda start, stop: ','.join(str(m).split(',')[start:stop])
        first_slice = slice_m(0,len(fit))+','+slice_m(len(fit), len(fit)+1).split(':')[0]
        second_slice = ','.join(m.split(',')[len(fit):])[2:]
        fun_args = [a.strip() for a in str(m).split(',')[1:len(fit)+1]]
        fun_args[-1] = fun_args[-1][0]
        for a,b in zip(fun_args,fit):
            second_slice = second_slice.replace(a.strip(),str(round(b,5)))
        self.legend += [prefix+ str(m) + 5*' ' + second_slice + suffix]


class Manipulator():
    def __init__(self, df):
        self._df = df

    ## Building block of Study class
    # data manipulation tools
  
    def assert_structure(self, df, structure):
        assert structure in df.columns, f"Structure {structure} not found"

    def define_index(self, df, samp, construct, cluster, index):
        if index in ['all','full'] :
            return df.index
        if index == 'roi':
            assert [roi in df.columns for roi in ['ROI_start','ROI_stop']], 'ROI_start and ROI_stop not found'
            return list(range(int(self.get_series(self._df, samp, construct, cluster)['ROI_start']), int(self.get_series(self._df, samp, construct, cluster)['ROI_stop'])))
        if type(index) in [list,tuple]:
            assert [i in list(range(len(self.get_series(self._df, samp, construct, cluster)['sequence']))) for i in index], 'Index out of range'
            return index
        if type(index) == str and not sum([int(a not in ['A','C','G','T']) for a in index]):
            return self.find_sequence_index(samp, construct, index)
        raise ValueError(f"Index {index} not recognized")

    def filter_base_paired(self, df_loc, base_paired):
        if base_paired != None:
            df_loc = df_loc[df_loc['paired'] == base_paired]
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

    def filter_flank(self, df, flank):
        if flank == None:
            return df
        if type(flank) == str:
            assert len((df_out:= df[df['flank'].apply(lambda x: flank == x)]).index) >=1, f"No row found with flank {flank}"
            return df_out
        if type(flank) in [list, tuple]:
            assert len((df_out:= df[df['flank'].apply(lambda x: x in flank)]).index) >=1, f"No row found with flanks {flank}"
            return df_out
        raise ValueError(f"Flank {flank} not recognized")        

    def filter_sub_lib(self, df, sub_lib):
        if sub_lib == None:
            return df
        if type(sub_lib) == str:
            assert len((df_out:= df[df['sub-library'].apply(lambda x: sub_lib == x)]).index) >=1, f"No row found with sub-library {sub_lib}"
            return df_out
        if type(sub_lib) in [list, tuple]:
            assert len((df_out:= df[df['sub-library'].apply(lambda x: x in sub_lib)]).index) >=1, f"No row found with sub-libraries {sub_lib}"
            return df_out
        raise ValueError(f"Sub-library {sub_lib} not recognized")

    def assert_SCC_exists(self, samp, construct, cluster):
        if samp != None:
            assert samp in list(self._df['samp']), f"Sample {samp} not found"
            assert construct in list(self._df[self._df['samp']==samp].construct.unique()), f"Construct {construct} not found"
            if cluster != None:
                assert len(self._df[(self._df['samp'] == samp) & (self._df['construct'] == construct) & (self._df['cluster'] == cluster)]) >= 1, f"No row found for {samp} {construct} {cluster}"
        else:
            assert construct in list(self._df.construct.unique()), f"Construct {construct} not found"
            if cluster != None:
                assert len(self._df[(self._df['construct'] == construct) & (self._df['cluster'] == cluster)]) >= 1, f"No row found for {samp} {construct} {cluster}"
            

    def find_sequence_index(self, samp, construct, sequence):
        df = self._df
        self.assert_SCC_exists(samp, construct, cluster = None)
        full_seq = df[(df['samp']==samp)&(df['construct']==construct)]['sequence'].iloc[0]
        assert (count:=full_seq.count(sequence)) == 1, f"{count} sequences {sequence} found in sequence of sample {samp} construct {construct} (sequence: {full_seq})"
        ind = full_seq.find(sequence)
        return list(range(ind, ind+len(sequence)))


    def get_col_across_constructs(self, samp:str, col:str, constructs='all', cluster=None, structure=None, base_type = ['A','C','G','T'], index='all', base_paired=None, flank=None, sub_lib=None ):
        args = locals()
        df_dont_touch = self._df
        df = self.filter_flank(df_dont_touch, flank)
        df = self.filter_sub_lib(df, sub_lib)

        stack = pd.DataFrame()
        if constructs == 'all':
            constructs = list(df[df.samp == samp].construct.unique())
        
        for c in constructs:
            stack = pd.concat((stack, self.get_SCC(construct=c, cols=[col], **{k:v for k,v in args.items() if k in self.get_SCC.__code__.co_varnames and k not in ['self','col']}).T))
        stack.index = constructs
        return stack


    def get_cols_by_sequence_across_constructs(self, samp:str, col:str, constructs='all', cluster=None, structure=None, base_type = ['A','C','G','T'], index='all', base_paired=None, flank=None, sub_lib=None ):
        args = locals()
        df_dont_touch = self._df
        df = self.filter_flank(df_dont_touch, flank)
        df = self.filter_sub_lib(df, sub_lib)
        assert (type(index) == str and not sum([int(a not in ['A','C','G','T']) for a in index])), f"Index is {index}, should be a str of 'A','C','G','T'."
        stack = pd.DataFrame()
        if constructs == 'all':
            constructs = list(df[df.samp == samp].construct.unique())
        
        for c in constructs:
            temp = self.get_SCC(construct=c, cols=[col], **{k:v for k,v in args.items() if k in self.get_SCC.__code__.co_varnames and k not in ['self','col']}).T
            temp.columns = list(index)
            stack = pd.concat((stack, temp))
        stack.index = constructs
        return stack




    def get_SCC(self, samp, construct, cols, cluster=0, structure=None, base_type = ['A','C','G','T'], index='all', base_paired=None,can_be_empty=False):
        """Returns a dataframe containing the content of a cluster of a sample-construct.

        Args:
            df (pd.Dataframe): A study dataframe.
            samp (str): The sample name.
            construct (str): The construct name.
            cols (list): The columns to be returned.
            cluster (int, optional): The cluster number. Defaults to 0.
            structure (str, optional): Structure to use for the 'paired' column, such as 'structure_ROI_DMS'. Defaults to 'structure'.
            base_type (list, optional): Bases to include. Defaults to ['A','C','G','T'].
            index (str, optional): Index to include. Defaults to 'all'. Can be a series of 0-indexes (ex: [43,44,45,48]), 'roi', 'all', or a unique sequence (ex: 'ATTAC')
            base_paired (_type_, optional): Base pairing to include. None is paired + unpaired, True is paired, False is unpaired. Defaults to None.
            can_be_empty (bool, optional): If True, returns an empty dataframe if no row is found. Defaults to False.

        Returns:
            _type_: dataframe containing the content of a cluster of a sample-construct.
        """

        self.assert_SCC_exists( samp, construct, cluster)
        df = self._df

        if structure != None:
            self.assert_structure(df, structure)
            assert len([c for c in cols if c.startswith('structure')])==0, f"{structure} should be entered as a structure= argument, not as a col. Cols can't contain structure."
            cols += [structure]

        remove_bases_flag = False
        if 'sequence' not in cols:
            cols += ['sequence']
            remove_bases_flag = True

        for col in cols:
            assert col in df.columns, f"Column {col} not found"

        df_loc = self.get_series(df, samp, construct, cluster, can_be_empty)
        for col in [c for c in cols if type(df_loc[c]) in [str]]:
            df_loc[col] = list(df_loc[col])

        df_loc = pd.DataFrame({col: df_loc[col] for col in cols})
        for st in [col for col in cols if col.startswith('structure')]:
            df_loc['paired'] = [{'.':False,'(':True,')':True}[x] for x in df_loc[st]]
            df_loc = df_loc.drop(columns=st)
        
        df_loc = df_loc.rename(columns={'sequence':'base'})
        index = self.define_index(df_loc, samp, construct, cluster, index)
        df_loc = self.filter_bases(df_loc, base_type, index, base_paired)
        if remove_bases_flag:
            df_loc = df_loc.drop(columns='base')
        return df_loc.sort_index()

    def get_series(self, df, samp, construct, cluster, can_be_empty=False):
        if samp == None:
            if cluster == None:
                df_out = df[df['construct'] == construct]
            else:
                df_out = df[(df['construct'] == construct) & (df['cluster'] == cluster)]            
        else:
            if cluster == None:
                df_out = df[(df['samp'] == samp) & (df['construct'] == construct)]
            else:
                df_out = df[(df['samp'] == samp) & (df['construct'] == construct) & (df['cluster'] == cluster)]
            assert len(df_out) <= 1, 'More than one row found'
        if not can_be_empty:
            assert len(df_out) >= 1, 'No row found'
        else: 
            if len(df_out) == 0: 
                return None
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
      

    def collect_x_y_paired_unpaired(self, cols:Tuple, samp, structure, cluster=None, max_mutation=MAX_MUTATION, base_type=['A','C','G','T'], index='all'):

        args = locals()
        args.pop('self')
        stack = {True:{'x':[],'y':[]},False:{'x':[],'y':[]}}
        df = self._df.copy()

        def stack_up(x,y,stack,is_paired):
            if len(x) != 0:
                stack[is_paired]['x'] += list(x)
                stack[is_paired]['y'] += list(y)

        def clean_stack(stack, max_mutation):
            for is_paired in [True,False]:
                for ax in ['x','y']:
                    stack[is_paired][ax] = np.array(stack[is_paired][ax]).reshape(1,-1)[0]
                if max_mutation:
                    assert 'mut_rates' in cols, 'max_mutation requires mut_rates'
                    stack[is_paired]['y'][stack[is_paired]['y'] >= max_mutation] = np.nan

            for is_paired in [True,False]:
                for ax in ['x','y']:
                    v = stack[is_paired]['x']
                    stack[is_paired]['x'] = stack[is_paired]['x'][~np.isnan(stack[is_paired][ax])]
                    stack[is_paired]['y'] = stack[is_paired]['y'][~np.isnan(stack[is_paired][ax])]
            return stack

        for construct in df[df['samp']==samp]['construct'].unique():
            for cluster in df[df['construct'] == construct].cluster.unique() if cluster is None else [cluster]:
                df_SCC = self.get_SCC(cols=cols.copy(),construct=construct,can_be_empty=True, **{k:v for k,v in args.items() if k in self.get_SCC.__code__.co_varnames and (k not in ['cols','df'])})
                for is_paired in [True,False]:
                    stack_up(df_SCC[df_SCC.paired == is_paired][cols[0]], df_SCC[df_SCC.paired == is_paired][cols[1]], stack, is_paired)
        stack = clean_stack(stack, max_mutation)
        return stack

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
