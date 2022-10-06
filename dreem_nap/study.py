from genericpath import exists
from random import sample
from typing import List
from dreem_nap import plotter, manipulator, util
from dreem_nap.loader import df_from_local_files
import pandas as pd
from dreem_nap import deltaG


class Study(object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'])
    """

    attr_list = ['name','samples']

    def __init__(self, name:str=None, samples:List[str]=None) -> None:
        """Creates a Study object.

        Args:
            name (str, optional): Short description (<~20 char) of your study. Defaults to None.
            samples (List[str], optional): names of your study's samples. Defaults to None.

        Raises:
            f: if conditions don't match the samples, or if the unit is missing.

        Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'arbitrary unit', 'Just an example study')
        >>> study.description
        'Just an example study'
        """

        self.name = name
        self.samples = samples
        self.constructs = None
        self._df = None

    def get_df(self, **kwargs):
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
        return self._df

    def set_df(self, df):
        self._df = df
        return df
    
    def get_constructs(self, samp:str):
        return self._df[self._df['samp'] == samp]['construct'].unique()

    def get_genes(self, samp:str, construct:str):
        return self._df[(self._df['samp'] == samp) & (self._df['construct'] == construct)]['gene'].unique()

    def get_clusters(self, samp:str, construct:str, gene:str):
        return self._df[(self._df['samp'] == samp) & (self._df['construct'] == construct)& (self._df['gene'] == gene)]['cluster'].unique()
    

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
        return cls(di['name'], di['samples'])

       
    def load_studies(studies_file_path:str):
        return load_studies(studies_file_path)


    def load_df_from_local_files(self, path_to_data:str, min_cov_bases:int, filter_by='study', index='all', base_type = ['A','C','G','T'], base_paired=None, structure=None)->pd.DataFrame:
        sub_df = util.SubDF.from_locals(locals())
        df = self.set_df(df_from_local_files(path_to_data, min_cov_bases, self.samples, self.name, filter_by, sub_df))
        self.constructs = df['construct'].unique()
        return df

    def get_col_across_constructs(self, **kwargs )->pd.DataFrame:
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
        return manipulator.Manipulator(self._df).get_col_across_constructs(**kwargs)
    
    # Plots

    def deltaG_per_sample(self, **kwargs)->util.OutputPlot:
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            samp (str): Sample of your sample-construct-cluster.
            deltaG (str): DeltaG to use as x axis.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to 'structure'.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            flank (str, optional): Flank or list of flanks to filter by. Defaults to None.
            sub_lib (str, optional): Sub-library or list of sub-libraries to filter by. Defaults to None.
            max_mutation (float, optional): Maximum mutation rate to plot. Defaults to 0.15.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            savefile (str, optional): Path to save the plot. Defaults to None.

        Returns:
            OutputPlot: Figure and data of the output plot.
        """
        return plotter.deltaG_per_sample(self._df, **kwargs)

    
    def deltaG_per_base(self, **kwargs)->util.OutputPlot:
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            construct (str): Construct of your row.
            experimental_variable (str): x axis column value, must be a per-sample attribute.
            region (str): Region of your row.
            cluster (str): Cluster of your row.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to 'structure'.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            max_mutation (float, optional): Maximum mutation rate to plot. Defaults to 0.15.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            savefile (str, optional): Path to save the plot. Defaults to None.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """
        return plotter.deltaG_per_base(self._df, **kwargs)
        

    def mutation_histogram(self, **kwargs):
        """Plot the mutation rates as histograms.

        Args:
            samp (str): Sample of your row.
            construct (str): Construct of your row.
            region (str): Region of your row.
            cluster (int, optional): Cluster of your row. Defaults to 0. 
            show_ci (bool, optional): Show confidence interval on the histogram. Defaults to True.
            savefile (str, optional): Path to save the plot. Defaults to None.

        Raises:
            Exception: plot_type is not ``index`` or ``partition``.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """
        return plotter.mutation_histogram(self._df, **kwargs)

    def base_coverage(self, **kwargs):
        """Plot the base coverage of several constructs in a sample.

        Args:
            samp (str): Sample of your rows.
            constructs (List[str]): Constructs of your rows.
            region (str): Region of your row.
            cluster (int, optional): Cluster of your row. Defaults to 0. 
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``. Can be a series of 0-indexes (ex: [43,44,45,48]), 'roi', 'all', or a unique sequence (ex: 'ATTAC')
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            base_paired (bool, optional): Base-pairing predicition to plot. Defaults to None.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to 'structure'.
            show_ci (bool, optional): Show confidence interval on the histogram. Defaults to True.
            savefile (str, optional): Path to save the plot. Defaults to None.
    
        Returns:
            OutputPlot: Figure, axis and data of the output plot.

        """
        return plotter.base_coverage(self._df, **kwargs)


def load_studies(studies_file_path:str)->dict[str:Study]:
    """Read formatted file with samples, and turn it into a dataframe containing studies.

    Args:
        studies_file_path (str): path+title of the csv file containing the samples.

    Returns:
        (pd.DataFrame): studies of the csv file, indexed by study.
    """

    studies_dict, studies_data = {}, pd.read_csv(studies_file_path)

    for col in studies_data.groupby('name')[Study.attr_list]:
        solo_item = lambda x: x[0] if len(set(x)) == 1 else x  
        studies_dict[col[0]] = {attr: solo_item(list(col[1][attr])) for attr in (Study.attr_list)} 

    return {k:Study.from_dict(v) for k,v in studies_dict.items()}

