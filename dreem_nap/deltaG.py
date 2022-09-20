from dreem_nap.util import OutputPlot, MplAttr, SubDF
from dreem_nap import manipulator
import matplotlib.pyplot as plt
from itertools import cycle
from typing import Tuple, List

class DeltaG:
    def __init__(self, man, **kwargs)->None:
        self.__man = man

    def per_sample(self, sub_df:SubDF, mpl_attr:MplAttr, deltaG:str, flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], **kwargs)->OutputPlot:

        fit = manipulator.Fit() 
        assert deltaG in self.__man._df.columns, f"deltaG arg {deltaG} isn't in df columns"
        data = self.__man.collect_x_y_paired_unpaired(cols=[deltaG,'mut_rates'], sub_df=sub_df)

        out = OutputPlot(data=data, mpl_attr=mpl_attr)       

        for is_paired, color in zip([True,False],['g','r']):
            plt.plot(data[is_paired]['x'],data[is_paired]['y'], color+'.')

        for color, is_paired, prefix in zip(['g','r'],[True,False], ['Paired bases ','Unpaired bases ']):
            style = cycle(['-','--',':','-.'])
            for m in models:
                plt.plot(*fit.predict(data[is_paired]['x'],data[is_paired]['y'],  m, prefix=prefix), color=color, linestyle=next(style))
            
        plt.legend(['Paired bases data','Unpaired bases data'] + fit.legend )

        if len(sub_df.index) >50:
            sub_df.index = str(sub_df.index)[:50]+'... ]'
        plt.title('  '.join([f"{k}: {v}" for k,v in sub_df.__dict__.items() if \
            k not in ['base_type','self']\
            and v is not None]))
        plt.xlabel(f'DeltaG: {deltaG} [Kcal/mol]')
        plt.ylabel(f'Mutation rate')

        [getattr(plt, arg)(kwargs[arg]) for arg in kwargs if hasattr(plt, arg)] 
        return out
    

    def per_base(self, sub_df:SubDF, mpl_attr:MplAttr, deltaG:str, flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], **kwargs)->OutputPlot:

        fit = manipulator.Fit() 
        assert deltaG in self.__man._df.columns, f"deltaG arg {deltaG} isn't in df columns"
        