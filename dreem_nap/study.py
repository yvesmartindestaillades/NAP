from re import L
from typing import Tuple, List, Dict


class Study:
    def __init__(self, name:str=None, samples:List[str]=None, conditions:List[float]=None, conditions_unit:str= None, description:str = None) -> None:
        
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

    def to_dict(self):
        out_dict = {}
        for attr in self.attr_list:
            if getattr(self, attr) != None:
                out_dict[attr] = getattr(self, attr)
        return out_dict

    def from_dict(self, di):
        for attr in  self.attr_list:
            try: 
                di[attr]
            except: 
                di[attr]=None
        self.__init__(di['name'], di['samples'], di['conditions'], di['conditions_unit'], di['description'])
        return self
