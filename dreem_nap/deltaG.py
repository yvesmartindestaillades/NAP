from dreem_nap.util import OutputPlot, MplAttr, SubDF
from dreem_nap import manipulator, util
import matplotlib.pyplot as plt
from itertools import cycle
from typing import Tuple, List
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from dreem_nap.manipulator import Fit

def per_sample(df:pd.DataFrame, samp:str, deltaG:str='deltaG_min', structure:str='structure', index='all', base_type=['A','C','G','T'], flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], **kwargs)->OutputPlot:

    fit = manipulator.Fit()
    man = manipulator.Manipulator(df)
    assert deltaG in man._df.columns, f"deltaG arg {deltaG} isn't in df columns"
    data = pd.DataFrame()
    sub_df= SubDF.from_locals(locals())
    for row in df[df.samp==samp].itertuples():
        sub_df.update(construct = row.construct, cluster=row.cluster)
        data = pd.concat((data,man.get_SCC(cols=[deltaG,'mut_rates','construct'],sub_df=sub_df,can_be_empty=True).reset_index()))

    hover_attr = ['construct','index','mut_rates',deltaG]
    tra = {}
    for is_paired, prefix in zip([True,False], ['Paired ','Unpaired ']):
        tra[prefix] = go.Scatter(
            x=data[data.paired == is_paired][deltaG], 
            y=data['mut_rates'], 
            text = data[hover_attr],
            mode='markers',
            name=prefix,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)])
        )         
        for m in models:
            x = data[(data.paired==is_paired)][deltaG]
            y = data[(data.paired==is_paired)]['mut_rates']
            fit = Fit()
            x_sorted, pred_y_sorted = fit.predict(x,y,m, prefix)
            tra[fit.get_legend()] = go.Scatter(
                x=x_sorted,
                y=pred_y_sorted,
                mode='lines',
                name=fit.get_legend())

    layout = dict(title = 'Mutation rates of paired / unpaired residues vs the predicted energy of the molecule',
            xaxis= dict(title= 'DeltaG',ticklen= 5,zeroline= False),
            yaxis= dict(title= 'Mutation rate ',ticklen= 5,zeroline= False),
            )
    fig = dict(data = list(tra.values()), layout = layout)
    iplot(fig)
    plot(fig, filename = 'filename.html', auto_open=False)
    return OutputPlot(data, fig)
    

def per_base(df:pd.DataFrame, construct:str, experimental_variable:str, structure:str='structure', index='all', base_type=['A','C','G','T'], max_mutation:float= 0.15, models:List[str]=[], **kwargs)->OutputPlot:

    fit = manipulator.Fit()
    man = manipulator.Manipulator(df)
    data = pd.DataFrame()
    sub_df= SubDF.from_locals(locals())
    assert construct in list(df.construct), f"{construct=} isn't in the dataframe"
    assert experimental_variable in df.columns, f"{experimental_variable=} isn't in the study columns"
    hover_attr = ['base','mut_rates','samp',experimental_variable]
    for row in df[df.construct==construct].itertuples():
        sub_df.update(samp = row.samp, construct = row.construct, cluster=row.cluster)
        data = pd.concat((data,man.get_SCC(cols=hover_attr.copy(),sub_df=sub_df,can_be_empty=True).reset_index()))
    data = data.reset_index().rename(columns={'level_0':'index_subsequence'})
    hover_attr+=['paired']
    data['Marker'] = data['paired'].apply(lambda x: {True: 0, False:1}[x])
    tra = {}
    for idx, row in data.groupby('index_subsequence'):
        name = f"({row['base'].iloc[0]},{row['index'].iloc[0]})"
        tra[idx] = go.Scatter(
            x= row['DMS_conc_mM'], 
            y= row['mut_rates'], 
            text = data[hover_attr],
            mode='markers',
            marker = dict(symbol = list(map(util.Setshape, data['Marker']))),
            name= name,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)])
        )         
        for m in models:
            x= row['DMS_conc_mM']
            y= row['mut_rates']
            fit = Fit()
            x_sorted, pred_y_sorted = fit.predict(x,y,m, name)
            
            tra[fit.get_legend()] = go.Scatter(
                x=x_sorted,
                y=pred_y_sorted,
                mode='lines',
                name=fit.get_legend())

    layout = dict(title = 'Mutation rates of paired / unpaired residues vs '+experimental_variable,
            xaxis= dict(title= experimental_variable,ticklen= 5,zeroline= False),
            yaxis= dict(title= 'Mutation rate ',ticklen= 5,zeroline= False),
            )

    fig = dict(data = list(tra.values()), layout = layout)
    iplot(fig)
    plot(fig, filename = 'filename.html', auto_open=False)
    return OutputPlot(data, fig)