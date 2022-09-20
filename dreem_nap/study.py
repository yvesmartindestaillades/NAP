from genericpath import exists
from random import sample
from typing import List
from dreem_nap import plotter, manipulator, util
from dreem_nap.loader import df_from_local_files
import pandas as pd


class Study(object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        conditions (List[float], optional): Values of the experimental condition that changes between the samples. Defaults to None.
        label (str, optional): Short description of the conditions. Defaults to None.
        description (str, optional): More information about your study. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'Example values [no unit]', 'Just an example study')
        >>> study.description
        'Just an example study'
        >>> study.to_dict()
        {'name': 'example', 'description': 'Just an example study', 'samples': ['A1', 'B2', 'B3'], 'label': 'Example values [no unit]', 'conditions': [10, 20, 30]}
        >>> di = {'name':'temperature','samples':['A1','B1','C3']}
        >>> study = Study.from_dict(di)
        >>> print(study.name, study.samples, study.description)
        temperature ['A1', 'B1', 'C3'] None
    """

    attr_list = ['name','description','samples','label','conditions']

    def __init__(self, name:str=None, samples:List[str]=None, conditions:List[float]=None, label:str= None, description:str = None) -> None:
        """Creates a Study object.

        Args:
            name (str, optional): Short description (<~20 char) of your study. Defaults to None.
            samples (List[str], optional): names of your study's samples. Defaults to None.
            conditions (List[float], optional): samples.csv column or values of the experimental condition that changes between the samples. Defaults to None.
            label (str, optional): Short description of the conditions. Defaults to None.
            description (str, optional): More information about your study. Defaults to None.

        Raises:
            f: if conditions don't match the samples, or if the unit is missing.

        Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'arbitrary unit', 'Just an example study')
        >>> study.description
        'Just an example study'
        """

        self.name = name
        self.description = description
        self.samples = samples
        self.constructs = None
        self.label = label
        self.conditions = conditions
        self._df = None
        self.plot = plotter.Plotter(self._df)
        self.mani = manipulator.Manipulator(self._df)

        if conditions != None and len(conditions) != len(samples):
            raise f"Number of samples ({len(samples)})and number of conditions ({len(conditions)}) don't match"

    def get_df(self, samp=None, construct=None, cluster=None, cols='all',structure=None,base_type=['A','C','G','T'],index='all',base_paired=None,sub_lib=None,flank=None,can_be_empty=False):
        """_summary_

        Returns:
            samp (str): The sample name.
            construct (str): The construct name.
            cols (list): The columns to be returned. Default is 'all'
            cluster (int, optional): The cluster number. Defaults to 0.
            structure (str, optional): Structure to use for the 'paired' column, such as 'structure_ROI_DMS'. Defaults to 'structure'.
            base_type (list, optional): Bases to include. Defaults to ['A','C','G','T'].
            index (str, optional): Index to include. Defaults to 'all'. Can be a series of 0-indexes (ex: [43,44,45,48]), 'roi', 'all', or a unique sequence (ex: 'ATTAC')
            base_paired (bool, optional): Base pairing to include. None is paired + unpaired, True is paired, False is unpaired. Defaults to None.
            sub_lib (str, optional): Sub-library
            flank (str, optional): Flank
            can_be_empty (bool, optional): If True, returns an empty dataframe if no row is found. Defaults to False.
        """
        df = self._df.copy()
        sub_df = util.SubDF.from_locals(locals())
        df = self.mani.filter_flank(df, flank)
        df = self.mani.filter_sub_lib(df, sub_lib)
        stack_index = False
        if samp == None:
            samp = df.samp.unique()
        if type(samp) in [str,int]:
            samp = [samp]
            stack_index = []
        if construct == None:
            construct = df.construct.unique()
        if type(construct) in [str, int]:
            construct = [construct]
        filter_by_cluster = cluster != None
        stack = pd.DataFrame()
        if cols=='all':
            cols = df.columns

        def stack_up(stack, sub_df, stack_index, cols):
            if stack_index:
                stack_index += [s+' - '+c]
            temp = self.mani.get_SCC(cols=cols,sub_df=sub_df,can_be_empty=True)
            if not temp.empty:
                stack = pd.concat((stack, pd.DataFrame({k:list(temp[k]) for k in cols})))
            return stack

        for s in list(set(samp)&set(df.samp.unique())) :
            for c in list(set(construct)&set(df.construct.unique())):
                sub_df.samp = s
                sub_df.construct = c
                if filter_by_cluster:
                    for cl in list(set(cluster)&set(df[(df.samp==s) & (df.construct==c)]['cluster'].unique())):
                        sub_df_temp = sub_df
                        sub_df_temp.cluster = cl
                        stack = stack_up(stack, sub_df_temp, stack_index, cols)
                else:
                    stack = stack_up(stack, sub_df, stack_index, cols)
        if stack_index:
            stack.index = stack_index
        return stack

    def set_df(self, df):
        self._df = df
        self.plot = plotter.Plotter(df)
        self.mani = manipulator.Manipulator(df)
        return df
    
    def get_constructs(self, samp:str):
        return self._df[self._df['samp'] == samp]['construct'].unique()

    def get_clusters(self, samp:str, construct:str):
        return self._df[(self._df['samp'] == samp) & (self._df['construct'] == construct)]['cluster'].unique()
        
    @classmethod
    def from_dict(cls, di:dict):
        """Set attributes of this Study object from a dictionary.

        Args:
            di (dict): a dictionary containing keys such as ['name','description','samples','label','conditions'].

        Returns:
            Study: a study object.

        Example:
        >>> di = {'name':'temperature','samples':['A1','B2','B3']}
        >>> study = Study().from_dict(di)
        >>> print(study.name, study.samples)
        temperature ['A1', 'B2', 'B3']
        """
        for attr in cls.attr_list:
            try: 
                di[attr]
            except: 
                di[attr]=None
        return cls(di['name'], di['samples'], di['conditions'], di['label'], di['description'])

       
    def load_studies(studies_file_path:str):
        return load_studies(studies_file_path)


    def load_df_from_local_files(self, path_to_data:str, min_cov_bases:int, filter_by='study', index='all', base_type = ['A','C','G','T'], base_paired=None, structure=None)->pd.DataFrame:
        sub_df = util.SubDF.from_locals(locals())
        df = self.set_df(df_from_local_files(path_to_data, min_cov_bases, self.samples, self.name, filter_by, sub_df))
        self.constructs = df['construct'].unique()
        return df

    def get_col_across_constructs(self, samp:str, col:str, construct='all', cluster=None, structure=None, base_type = ['A','C','G','T'], index='all', base_paired=None, flank=None, sub_lib=None )->pd.DataFrame:
        """Returns a dataframe containing the column col for provided constructs in a sample

        Args:
            samp (str): Sample(s) of your sample-construct-cluster. A single sample or a list of samples.
            col (list): The column to be returned.
            constructs (str): The constructs name. Defaults to 'all'.
            cluster (int, optional): The cluster number. Defaults to 0.
            structure (str, optional): Structure to use for the 'paired' column, such as 'structure_ROI_DMS'. Defaults to 'structure'.
            base_type (list, optional): Bases to include. Defaults to ['A','C','G','T'].
            index (str, optional): Index to include. Defaults to 'all'. Can be a series of 0-indexes (ex: [43,44,45,48]), 'roi', 'all', or a unique sequence (ex: 'ATTAC')
            base_paired (bool, optional): Base pairing to include. None is paired + unpaired, True is paired, False is unpaired. Defaults to None.
            flank (str or list, optional): Flank or list of flanks to filter constructs by. Defaults to None.
            sub_lib (str or list, optional): Sub-library or list of sub-libraries to filter constructs by. Defaults to None.

        Returns:
            pd.Dataframe: content of the column col across constructs. Columns names are the indexes provided by index and index names (y axis) are the constructs.
        """
        sub_df = util.SubDF.from_locals(locals)
        return self.__mani.get_col_across_constructs(sub_df, col, flank, sub_lib )


def load_studies(studies_file_path:str)->dict[str:Study]:
    """Read formatted file with samples, and turn it into a dataframe containing studies.

    Args:
        studies_file_path (str): path+title of the csv file containing the samples.

    Returns:
        (pd.DataFrame): studies of the csv file, indexed by study.

    Example:
        >>> from dreem_nap import data_wrangler     
        >>> study_file = 'samples.csv'
        >>> samples = ['"name","description","samples","conditions","label"',
        ... ',,,,',
        ... '"salt","Change the Na concentration","A6",0.15,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","B6",0.3,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","C6",0.6,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","D6",1,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","E6",1.2,"Na quantity [M]"',
        ... ',,,,',
        ... '"spermidine","Change the Spermidine concentration","B3",0.01,"Spermidine quantity [mM]"',
        ... '"spermidine","Change the Spermidine concentration","D3",1,"Spermidine quantity [mM]"',
        ... '"spermidine","Change the Spermidine concentration","E3",10,"Spermidine quantity [mM]"',
        ... '"spermidine","Change the Spermidine concentration","F3",100,"Spermidine quantity [mM]"']
        >>> with open(study_file, 'w') as f:
        ...     f.writelines(f"{l}\\n" for l in samples)
        >>> # Your studies under the shape of a dataframe
        >>> df_studies = data_wrangler.load_studies(study_file) 
        >>> df_studies
                                name                           description                     samples                          label                        conditions
        salt                    salt           Change the Na concentration        [A6, B6, C6, D6, E6]                Na quantity [M]        [0.15, 0.3, 0.6, 1.0, 1.2]
        spermidine        spermidine   Change the Spermidine concentration            [B3, D3, E3, F3]       Spermidine quantity [mM]          [0.01, 1.0, 10.0, 100.0]
        >>> temp = df_studies.to_dict(orient='index')
        >>> # Your studies under the shape of a dictionary of Study
        >>> studies = {study: Study.from_dict(temp[study])  for study in temp}
        >>> print(f"Here are the studies: {studies.keys()}")
        dict_keys(['salt', 'spermidine'])
        >>> study_name = 'salt' 
        >>> study = studies[study_name] 
        >>> print(f"Here is your study {study.to_dict()}" )
        Here is your study {'name': 'salt', 'description': 'Change the Na concentration', 'samples': ['A6', 'B6', 'C6', 'D6', 'E6'], 'label': 'Na quantity [M]', 'conditions': [0.15, 0.3, 0.6, 1.0, 1.2]}
    """

    studies_dict, studies_data = {}, pd.read_csv(studies_file_path)

    for col in studies_data.groupby('name')[Study.attr_list]:
        solo_item = lambda x: x[0] if len(set(x)) == 1 else x  
        studies_dict[col[0]] = {attr: solo_item(list(col[1][attr])) for attr in (Study.attr_list)} 

    return {k:Study.from_dict(v) for k,v in studies_dict.items()}

