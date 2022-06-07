from typing import Tuple, List, Dict

class study:
    def __init__(self, name:str, samples:List[str], conditions:List[float]=None, condition_unit:str= None, description:str = None) -> None:
        self.name = name
        self.description = description
        self.samples = samples
        self.condition_unit = condition_unit
        self.conditions = conditions

        if conditions != None and len(conditions) != len(samples):
            raise f"Number of samples ({len(samples)})and number of conditions ({len(conditions)}) don't match"
        
        if conditions != None and condition_unit == None:
            raise "Conditions were given without a unit."

    def to_dict(self):
        out_dict = {}
        for attr in ['name','description','samples','condition_unit','conditions']:
            if getattr(self, attr) != None:
                out_dict[attr] = getattr(self, attr)
        return out_dict