from binascii import a2b_hex
from typing import Tuple, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath(""))
from dreem_nap.manipulator import Manipulator

from dreem_nap.clustering import Clustering
from dreem_nap.util import *
from dreem_nap import deltaG
import plotly

from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from dreem_nap.manipulator import Fit, Manipulator

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


def mutation_histogram(df, samp:str, construct:str, region:str=None, cluster:int=None, structure:str=None, show_ci:bool=True, savefile=None)->OutputPlot:

    mh = Manipulator(df).get_series(df, SubDF.from_locals(locals()))
    xaxis_coordinates = [i for i in range(len(mh.sequence) -1)]

    mut_y = []
    for pos in range(len(mh.sequence)):
        try:
            mut_frac = mh.mut_bases[pos] / mh.info_bases[pos]
        except:
            mut_frac = 0.0
        mut_y.append(mut_frac)
        mut_frac = round(mut_frac, 5)

    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    colors = []
    ref_bases = []
    hover_attr = pd.DataFrame({'mut_rate':mut_y,
                                'base':list(mh.sequence), 
                                'index':list(range(len(mh.sequence))),
                                'paired':[{'.':True, '(':False,')':False}[s] for s in mh.structure]})
    for i in range(len(mh.sequence)):
        if i >= len(mh.sequence)-1:
            continue
        colors.append(cmap[mh.sequence[i - 1]])
        ref_bases.append(mh.sequence[i - 1])
    mut_trace = go.Bar(
            x=xaxis_coordinates,
            y=mut_y,
            text=hover_attr,
            marker=dict(color=colors),
            showlegend=False,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)]),
        )   
    
    if show_ci:
        mut_trace.update(
           error_y=dict(
                type='data',
                symmetric=False,
                array=mh.poisson_high,
                arrayminus=mh.poisson_low
                )
        )

    mut_fig_layout = go.Layout(
            title=f"{mh.samp} - {mh.construct} - {mh.cluster}",
            xaxis=dict(title="Bases"),
            yaxis=dict(title="Mutation rate", range=[0, 0.1]),
            plot_bgcolor="white"

    )
    mut_fig = go.Figure(data=mut_trace, layout=mut_fig_layout)
    seqs = list(mh.sequence)
    if mh.structure is not None:
        db = list(mh.structure)
    else:
        db = " " * len(seqs)
    mut_fig.update_yaxes(
            gridcolor='lightgray',
            linewidth=1,
            linecolor='black',
            mirror=True
    )
    mut_fig.update_xaxes(
            linewidth=1,
            linecolor='black',
            mirror=True
    )
    mut_fig.update_xaxes(
            tickvals=xaxis_coordinates,
            ticktext=["%s<br>%s" % (x, y) for (x, y) in zip(seqs, db)],
            tickangle=0
    )
    iplot(mut_fig)

    if savefile != None:
        plot(mut_fig, filename = savefile, auto_open=False)

    return OutputPlot(mh, mut_fig)

def deltaG_per_sample(df:pd.DataFrame, samp:str, deltaG:str='deltaG_min', structure:str='structure', index='all', base_type=['A','C','G','T'], flank:str=None, sub_lib:str=None, max_mutation:float= 0.15, models:List[str]=[], savefile=None)->OutputPlot:

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
    if savefile != None:
        plot(fig, filename = savefile, auto_open=False)
    return OutputPlot(data, fig)
    
def deltaG_per_base(df:pd.DataFrame, construct:str, region:str, cluster:int, experimental_variable:str, structure:str='structure', index='all', base_type=['A','C','G','T'], max_mutation:float= 0.15, models:List[str]=[], savefile=None)->OutputPlot:

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
            mode='lines+markers',
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
    if savefile != None:
        plot(fig, filename = savefile, auto_open=False)
    return OutputPlot(data, fig)



def base_coverage(df, samp:str, constructs:str='all', gene:str=None, cluster:int=None, savefile='base_coverage.html')->OutputPlot:
    if constructs == 'all':
        constructs = list(df[df.samp==samp]['construct'].unique())
    trace = [
        go.Scatter(
            x= np.array([i for i in range(len(df.sequence.iloc[0]))]),
            y= np.array([df['min_cov_bases'].iloc[0] for i in range(len(df.sequence.iloc[0]))]) , 
            name='Min coverage bases',
            mode='lines')
    ]
    for construct in constructs:
        mh = Manipulator(df).get_series(df, SubDF.from_locals(locals()))
        x = np.array([i for i in range(len(mh.sequence))])
        y = np.array([int(mh.info_bases[i]) for i in range(len(mh.sequence))])
        trace.append(go.Scatter(
            x=x,
            y=y, 
            name=construct,
            mode='lines'))

    layout = go.Layout(
        title=f"{mh.samp} - {mh.construct} - {mh.cluster}",
        xaxis=dict(title="Bases"),
        yaxis=dict(title="Info bases"),
        plot_bgcolor="white"
        )

    fig = go.Figure(data=trace, layout=layout)
    iplot(fig)

    if savefile != None:
        plot(fig, filename = savefile, auto_open=False)

    return OutputPlot(data=pd.DataFrame({t['name']:{'x':t['x'], 'y':t['y']} for t in trace}).T, fig=fig)

