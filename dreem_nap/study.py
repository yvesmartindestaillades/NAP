from genericpath import exists
from random import sample
from typing import List
from dreem_nap import manipulator, util, plotter
from dreem_nap.loader import df_from_local_files
import pandas as pd
from dreem_nap import deltaG

NOT_IN_GET_DF = ['show_ci','savefile','deltaG','models','experimental_variable','auto_open']

class Study(object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'])
    """

    attr_list = ['name','samples']

    def __init__(self, path_to_data, samples=None, min_cov_bases=0, filter_by='sample') -> None:
        """Creates a Study object.

        Args:
            samples (List[str], optional): List of samples to load. Defaults to None.
            min_cov_bases (int, optional): Minimum number of base coverage for a row to be filtered-in. Defaults to 0.
            filter_by (str, optional): Filter rows by sample or study. When filtered by study, if a row passes the filter, rows with the same 'construct', 'section' and 'cluster' fields for all other samples have a sufficient base coverage. Defaults to 'sample'.            

        Example:
            >>> study = Study(path_to_data='data/my_study.csv', 
                              samples=['A1', 'B2', 'B3'], 
                              min_cov_bases=1000, 
                              filter_by='study')
        """

        self.samples = samples
        self.df = pd.read_csv(path_to_data)
        self.df = self.df[self.df['worst_cov_bases'] >= min_cov_bases]
        for col in [ 'mut_bases', 'info_bases','del_bases','ins_bases','cov_bases','mut_rates'] + \
            [c for c in self.df.columns.tolist() if (c.startswith('mod_bases') or c.startswith('poisson'))]:
            self.df[col] = self.df[col].apply(lambda x: [float(b) for b in x[1:-1].split(' ') if b != '' and b != '\n'])
        self.df = manipulator.get_df(df=self.df, sample=samples, min_cov_bases=min_cov_bases)
        if filter_by == 'study':
            self.df = manipulator.filter_by_study(self.df)

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

    def get_df(self, **kwargs):
        return manipulator.get_df(self.df, **kwargs)

    def get_samples(self):
        return self.df.sample.unique()

    def get_constructs(self, sample:str):
        return self.df[self.df['sample'] == sample]['construct'].unique()

    def get_genes(self, sample:str, construct:str):
        return self.df[(self.df['sample'] == sample) & (self._df['construct'] == construct)]['section'].unique()

    def get_clusters(self, sample:str, construct:str, section:str):
        return self.df[(self.df['sample'] == sample) & (self._df['construct'] == construct)& (self._df['section'] == section)]['cluster'].unique()
       
    def load_studies(studies_file_path:str):
        return load_studies(studies_file_path)


    def mutation_histogram(self, **kwargs):
        """Plot the mutation rates as histograms.
        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            construct (list, int, str, optional): Filter rows by construct (list of constructs or just a construct). Defaults to None.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
            base_index (list, int, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base type. Defaults to ['A','C','G','T'].
            base_pairing (bool, optional): Filter per-base attributes (mut_rates, sequence, etc) by predicted base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            show_ci(bool, optional): Show confidence intervals. Defaults to True.
            savefile(str, optional): Path to save the plot. Defaults to None.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """

        return plotter.mutation_histogram(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k not in NOT_IN_GET_DF}), **{k:v for k,v in kwargs.items() if k in plotter.mutation_histogram.__code__.co_varnames})

    def deltaG_per_sample(self, **kwargs)->util.OutputPlot:
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
            construct (list, int, str, optional): Filter rows by construct (list of constructs or just a construct). Defaults to None.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
            min_cov_bases (int, optional): Filter rows by a minimum threshold for base coverage. Defaults to 0.
            base_index (list, int, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base type. Defaults to ['A','C','G','T'].
            base_pairing (bool, optional): Filter per-base attributes (mut_rates, sequence, etc) by predicted base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            show_ci(bool, optional): Show confidence intervals. Defaults to True.
            savefile(str, optional): Path to save the plot. Defaults to None.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            **kwargs: Additional arguments to pass to filter rows by. Ex: flank='flank_1' will keep only rows with flank=flank_1. 

        Returns:
            OutputPlot: Figure and data of the output plot.
        """
        return plotter.deltaG_per_sample(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k not in NOT_IN_GET_DF}), **{k:v for k,v in kwargs.items() if k in plotter.deltaG_per_sample.__code__.co_varnames})

    
    def variable_exp_across_samples(self, **kwargs)->util.OutputPlot:
        """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.

        Args:
            construct (str): Construct of your row.
            experimental_variable (str): x axis column value, must be a per-sample attribute.
            section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
            cluster (str): Cluster of your row.
            structure (str, optional): Structure to use for base_paired filtering. Defaults to 'structure'.
            index (_type_, optional): Indexes to plot. Defaults to ``'all'``.
            base_type (List[str], optional): Bases type to plot. Defaults to ``['A','C','G','T']``.
            max_mutation (float, optional): Maximum mutation rate to plot. Defaults to 0.15.
            RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing and predicting deltaG. Defaults to False.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            savefile (str, optional): Path to save the plot. Defaults to None.

        Returns:
            OutputPlot: Figure, axis and data of the output plot.
        """
        return plotter.variable_exp_across_samples(manipulator.get_df(self.df, **{k:v for k,v in kwargs.items() if k not in NOT_IN_GET_DF}), **{k:v for k,v in kwargs.items() if k in plotter.variable_exp_across_samples.__code__.co_varnames})
        


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
        return 0# plotter.base_coverage(self._df, **kwargs)



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
