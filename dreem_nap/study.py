from re import L
from typing import Tuple, List, Dict
from dreem_nap import data, plot, data_manip
import pandas as pd

class Study(data.Data, plot.Plot, data_manip.Data_manip,  object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        conditions (List[float], optional): Values of the experimental condition that changes between the samples. Defaults to None.
        title (str, optional): Short description of the conditions. Defaults to None.
        description (str, optional): More information about your study. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'Example values [no unit]', 'Just an example study')
        >>> study.description
        'Just an example study'
        >>> study.to_dict()
        {'name': 'example', 'description': 'Just an example study', 'samples': ['A1', 'B2', 'B3'], 'title': 'Example values [no unit]', 'conditions': [10, 20, 30]}
        >>> di = {'name':'temperature','samples':['A1','B1','C3']}
        >>> study = Study().from_dict(di)
        >>> print(study.name, study.samples, study.description)
        temperature ['A1', 'B1', 'C3'] None
    """

    def __init__(self, name:str=None, samples:List[str]=None, conditions:List[float]=None, title:str= None, description:str = None) -> None:
        """Creates a Study object.

        Args:
            name (str, optional): Short description (<~20 char) of your study. Defaults to None.
            samples (List[str], optional): names of your study's samples. Defaults to None.
            conditions (List[float], optional): values of the experimental condition that changes between the samples. Defaults to None.
            title (str, optional): Short description of the conditions. Defaults to None.
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
        self.title = title
        self.conditions = conditions
        self.attr_list = ['name','description','samples','title','conditions']

        if conditions != None and len(conditions) != len(samples):
            raise f"Number of samples ({len(samples)})and number of conditions ({len(conditions)}) don't match"


    def to_dict(self)->dict:
        """Casts the Study object into a dictionary.

        Returns:
            dict: a dictionary form of the Study object.

        Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'Example values [no unit]', 'Just an example study')
        >>> study.to_dict()
        {'name': 'example', 'description': 'Just an example study', 'samples': ['A1', 'B2', 'B3'], 'title': 'Example values [no unit]', 'conditions': [10, 20, 30]}
         """

        out_dict = {}
        for attr in self.attr_list:
            if getattr(self, attr) != None:
                out_dict[attr] = getattr(self, attr)
        return out_dict

    def from_dict(self, di:dict):
        """Set attributes of this Study object from a dictionary.

        Args:
            di (dict): a dictionary containing keys such as ['name','description','samples','title','conditions'].

        Returns:
            Study: a study object.

        Example:
        >>> di = {'name':'temperature','samples':['A1','B2','B3']}
        >>> study = Study().from_dict(di)
        >>> print(study.name, study.samples)
        temperature ['A1', 'B2', 'B3']
        """
        for attr in  self.attr_list:
            try: 
                di[attr]
            except: 
                di[attr]=None
        self.__init__(di['name'], di['samples'], di['conditions'], di['title'], di['description'])
        return self


if __name__ == '__main__':
    temp = Study('temperature',['A1','B2','B3'], [10, 20, 30], 'Example values [no unit]', 'Just an example study')
    temp.create_df('data/DEMULTIPLEXED',10)
    print(temp.df)
    fig = temp.mut_histogram(temp.samples[0], '9572', 'index')
    