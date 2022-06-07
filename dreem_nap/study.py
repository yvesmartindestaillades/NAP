import re


class study:
    def __init__(self, name:str, samples:list(str), unit:str= None, values:list(float)=None, description:str = None) -> None:
        self.name = name
        self.description = description
        self.samples = samples
        self.unit = unit
        self.values = values

        if values != None and len(values) != len(samples):
            raise:  