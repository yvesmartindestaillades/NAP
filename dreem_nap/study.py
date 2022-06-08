from re import L
from typing import Tuple, List, Dict


class Study:
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        conditions (List[float], optional): Values of the experimental condition that changes between the samples. Defaults to None.
        conditions_unit (str, optional): Unit of the condition that changes between the samples. Defaults to None.
        description (str, optional): More information about your study. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'arbitrary unit', 'Just an example study')
        >>> study.description
        'Just an example study'
        >>> study.to_dict()
        {'name': 'example', 'description': 'Just an example study', 'samples': ['A1', 'B2', 'B3'], 'conditions_unit': 'arbitrary unit', 'conditions': [10, 20, 30]}
        >>> di = {'name':'temperature','samples':['A1','B1','C3']}
        >>> study = Study().from_dict(di)
        >>> print(study.name, study.samples, study.description)
        temperature ['A1', 'B1', 'C3'] None
    """

    def __init__(self, name:str=None, samples:List[str]=None, conditions:List[float]=None, conditions_unit:str= None, description:str = None) -> None:
        """Creates a Study object.

        Args:
            name (str, optional): Short description (<~20 char) of your study. Defaults to None.
            samples (List[str], optional): names of your study's samples. Defaults to None.
            conditions (List[float], optional): values of the experimental condition that changes between the samples. Defaults to None.
            conditions_unit (str, optional): unit of the condition that changes between the samples. Defaults to None.
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
        self.conditions_unit = conditions_unit
        self.conditions = conditions
        self.attr_list = ['name','description','samples','conditions_unit','conditions']

        if conditions != None and len(conditions) != len(samples):
            raise f"Number of samples ({len(samples)})and number of conditions ({len(conditions)}) don't match"
        
        if conditions != None and conditions_unit == None:
            raise "Conditions were given without a unit."

    def to_dict(self)->dict:
        """Casts the Study object into a dictionary.

        Returns:
            dict: a dictionary form of the Study object.

        Example:
        >>> study = Study('example',['A1', 'B2', 'B3'], [10, 20, 30], 'arbitrary unit', 'Just an example study')
        >>> study.to_dict()
        {'name': 'example', 'description': 'Just an example study', 'samples': ['A1', 'B2', 'B3'], 'conditions_unit': 'arbitrary unit', 'conditions': [10, 20, 30]}
        """

        out_dict = {}
        for attr in self.attr_list:
            if getattr(self, attr) != None:
                out_dict[attr] = getattr(self, attr)
        return out_dict

    def from_dict(self, di:dict):
        """Set attributes of this Study object from a dictionary.

        Args:
            di (dict): a dictionary containing keys such as ['name','description','samples','conditions_unit','conditions'].

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
        self.__init__(di['name'], di['samples'], di['conditions'], di['conditions_unit'], di['description'])
        return self
