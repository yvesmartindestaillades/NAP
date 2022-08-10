import pandas as pd
import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
from dreem_nap import utils



class Data_manip(object):
    def get_construct_attribute(self, column:str)->pd.DataFrame:   

        """Returns a column of the dataframe w.r.t constructs. 
        
        The value of this column must be only dependent on the construct, and not on the sample. 
        For example, the sequence can be read with this function, but not the mutation rate.  
        
        Args:
            attribute: the dataframe's column values that you want to look at.

        Returns:
            A dataframe with the constructs as index and the column as data.
        """   
        if not self.df.empty:
            return self.df.set_index('construct').sort_values(column)[column].groupby('construct').apply(lambda x:np.array(x)[0]).sort_values()
        


    def get_roi_info(self, samp:str, construct:int, bases_type:list[str]=['A','C'], structure = 'full', overlay = 0, roi_range=None)->pd.DataFrame:
        """Returns a dataframe of the ROI of a specific (samp, construct).

        Args:
            samp: a specific sample.
            construct: a specific construct.
            bases_type: list of the bases types to filter-in
            structure: 'full', 'roi', or 'both'. If 'full' or 'roi', the index 'paired' of the output will be corresponding to the structure prediction of the full RNA or only the ROI, respectively. If 'both', the output will be indexed w.r.t 'paired_full' and 'paired_roi'.  
            overlay (str or int or tuple[int]): extend/shrink roi
                'all': the roi is all bases
                int-type argument: the roi is the subsequence [start_roi_index-overlay, end_roi_index+overlay] 
                tuple[int]-type argument: the roi is the subsequence [start_roi_index-overlay[0], end_roi_index+overlay[1]] 
            roi_range: default is None. Array of base indexes (list[int]). ex: [80, 83, 91]. Base-0 index.

        Return:
            Indexes:
                base: A, C, G, T.
                paired: pairing prediction of an RNA structure prediction software.
                roi_structure_comparison: comparison between the pairing prediction of the ROI sequence alone and the entire sequence. 0 means no difference, 1 means difference. 
                index: base-0 index
            Columns:
                mut_rate: mutation rate of this base.
                roi_deltaG: the deltaG of the ROI sequence predicted b a RNA structure prediction software. 

        Raise:
            structure is 'roi' or 'full' and overlay is expanding the roi.
        """ 

        np.seterr(invalid='ignore')
        df = self.df
        df_SC = df.set_index(['samp','construct']).loc[(samp,construct)]

        assert not (overlay != 0 and roi_range != None), "overlay and roi_range are uncompatible arguments"

        roi_range, _ = utils.roi_range_calc(overlay, roi_range,                                                 
                                        roi_bounds=[df[df.construct==construct]['roi_start_index'].iloc[0], df[df.construct==construct]['roi_end_index'].iloc[0]],
                                        full_bounds=[df[df.construct==construct]['start'].iloc[0]-1, df[df.construct==construct]['end'].iloc[0]-1])
            
        assert not (structure != 'full' and (min(roi_range) < int(df_SC['roi_start_index']) or max(roi_range) > int(df_SC['roi_end_index']))), "Impossible to expand the roi when using roi-based structure prediction"

        df_roi = pd.DataFrame({'mut_rate':pd.Series(np.array(df_SC[f"mut_bases"][1:])/np.array(df_SC[f"info_bases"][1:]), dtype=object),
                                'base':list(df_SC['sequence']),
                                'paired': np.array([bool(x != '.') for x in list(df_SC['full_structure'])]),\
                                'base_pairing_prob': df_SC[f"base_pairing_prob"][1:],\
                                'roi_structure_comparison': pd.Series(list(df_SC['roi_structure_comparison']),index=list(range(df_SC['roi_start_index'], df_SC['roi_end_index'])))\
                                ,'roi_deltaG':df_SC['roi_deltaG']})\
                                .reset_index()
                                
        df_roi = df_roi.where(df_roi['base'].isin(bases_type))#.dropna()
        
        df_roi = df_roi[df_roi['index'].notnull()]
        df_roi['index'] =  df_roi['index'].astype(int)

        df_roi = df_roi[df_roi['index'].isin(roi_range)]

        if structure in ['roi','ROI','both']:
            df_roi['paired_roi'] = df_roi.apply(lambda row:  bool((int(row['paired'])+int(row['roi_structure_comparison']))%2)  , axis=1 )

        if structure == 'both':
            df_roi = df_roi.rename(columns={'paired':'paired_full'})

        if structure == 'roi':
            df_roi = df_roi.drop(columns=['paired'])
            df_roi = df_roi.rename(columns={'paired_roi':'paired'})

        if structure in ['roi','ROI','full']:
            df_roi = df_roi.set_index(['base', 'paired', 'roi_structure_comparison','index'])

        if structure == 'both':
            df_roi = df_roi.set_index(['base', 'paired_full', 'paired_roi', 'roi_structure_comparison','index'])

        return df_roi

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
        np.seterr(invalid='ignore')
        full_path = utils.make_path(path)
        df_print = self.df[self.df.samp.isin(samples)]
        df_print = df_print[columns] 
        np.set_printoptions(suppress=True)
        if 'mut_rate' in columns:
            df_print['mut_rate'] = df_print.apply(lambda row: np.float32(np.array(row['mut_bases'])/np.array(row['info_bases'])), axis=1)
        df_print = df_print.reset_index().drop(columns='index')
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
        all_samples = list(self.df.samp.unique())
        these_samples = np.array(all_samples)[np.random.randint(0, len(all_samples),n_samples)] 
        these_constructs = []
        for samp in these_samples:
            constructs = self.df[self.df.samp==samp].construct.unique()
            these_constructs.append(constructs[np.random.randint(0, len(constructs))])
        if n_samples == 1:
            these_samples = these_samples[0]
        if n_constructs == 1:
            these_constructs = these_constructs[0]

        return these_samples, these_constructs    

