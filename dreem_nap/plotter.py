from binascii import a2b_hex
from typing import Tuple, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath(""))

from dreem_nap.util import *
import plotly

from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

from dreem_nap.util import OutputPlot, MplAttr, SubDF
from dreem_nap import manipulator, util
import matplotlib.pyplot as plt
from itertools import cycle
from typing import Tuple, List
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import scipy
from scipy.optimize import curve_fit



class Fit(object):
    def __init__(self) -> None:
        self.legend = ''
    
    def get_legend(self):
        return self.legend

    def predict(self, x, y, model, prefix='', suffix=''):
        fit = self.fit(x,y,model)
        m = eval(model)
        try:
            linreg  = scipy.stats.linregress(y,m(x,*fit))
            self.rvalue = round(linreg.rvalue,5)
        except:
            self.rvalue = 'error'

        self._generate_legend(fit, model, prefix, suffix)
        return np.sort(x), m(np.sort(x),*fit)

    def fit(self, x,y, model):
        fit = curve_fit(eval(model), x, y)[0]
        return fit

    def _generate_legend(self, fit, m, prefix, suffix):
        slice_m = lambda start, stop: ','.join(str(m).split(',')[start:stop])
        first_slice = slice_m(0,len(fit))+','+slice_m(len(fit), len(fit)+1).split(':')[0]
        second_slice = ','.join(m.split(',')[len(fit):])[2:]
        fun_args = [a.strip() for a in str(m).split(',')[1:len(fit)+1]]
        fun_args[-1] = fun_args[-1][0]
        for a,b in zip(fun_args,fit):
            second_slice = second_slice.replace(a.strip(),str(round(b,5)))
        self.legend = prefix+ second_slice + suffix +f'\n R2={self.rvalue}'


def mutation_histogram(df, show_ci:bool=True, savefile=None, auto_open=False)->OutputPlot:
    assert len(df) == 1, "df must have only one row"
    mh = df.iloc[0]
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    
    traces, layouts = [], []
    mh_unrolled = pd.DataFrame({'mut_rate':list(mh.mut_rates), 'base':list(mh.sequence), 'filtered_index':list(mh.filtered_index), 'poisson_high':list(mh.poisson_high), 'poisson_low':list(mh.poisson_low), 'paired':list(mh.structure_selected)})

    for bt in set(mh['sequence']):
        df_loc = mh_unrolled[mh_unrolled['base'] == bt]
        if len(df_loc) == 0:
            continue

        hover_attr = pd.DataFrame({'mut_rate':list(df_loc.mut_rate),
                                        'base':list(df_loc.base), 
                                        'index':list(df_loc['filtered_index']),
                                        'paired':[{'.':True, '(':False,')':False}[s] for s in df_loc.paired]})
        traces.append(go.Bar(
            x= np.array(df_loc['filtered_index']),
            y= np.array(df_loc['mut_rate']),
            name=bt,
            marker_color=cmap[bt],
            text = hover_attr,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)]),
            ))
        if show_ci:
            traces[-1].update(
                        error_y=dict(
                        type='data',
                        symmetric=False,
                        array=df_loc['poisson_high'], 
                        arrayminus=df_loc['poisson_low']
                        ))

    
        mut_fig_layout = go.Layout(

        )

    fig = go.Figure(data=traces, layout=mut_fig_layout)

    fig.update_layout(title=f"{mh['sample']} - {mh['construct']} - {mh['section']} - {mh['cluster']}",
                        xaxis=dict(title="Sequence"),
                        yaxis=dict(title="Mutation rate", range=[0, 0.1]))
   
    fig.update_yaxes(
            gridcolor='lightgray',
            linewidth=1,
            linecolor='black',
            mirror=True,
            autorange=True
    )
    fig.update_xaxes(
            linewidth=1,
            linecolor='black',
            mirror=True,
            autorange=True
    )

    fig.update_xaxes(
            tickvals=mh_unrolled['filtered_index'],
            ticktext=["%s %s" % ({'.':'(P)','(':'(U)',')':'(U)'}[x], str(y)) for (x,y) in zip(mh['structure_selected'],mh['filtered_index'])],
            tickangle=90,
            autorange=True
    )

    iplot(fig)
    if savefile != None:
        plot(fig, filename = savefile, auto_open=auto_open)

    return OutputPlot(mh, fig)


def deltaG_per_sample(df:pd.DataFrame, models:List[str]=[], savefile=None, auto_open=False)->OutputPlot:

    df_temp = pd.DataFrame()
    for _, row in df.iterrows():
        df_temp = pd.concat([df_temp, pd.DataFrame({'construct':row.construct, 'index':row.filtered_index, 'mut_rates':row.mut_rates, 'deltaG':row['deltaG_selected'], 'paired':[s !='.' for s in row.structure_selected]}, index=list(range(len(row.filtered_index))))])
    df = df_temp.reset_index()

    hover_attr = ['construct','index','mut_rates','deltaG']
    tra = {}
    for is_paired, prefix in zip([True,False], ['Paired ','Unpaired ']):
        x=np.array(df[df.paired == is_paired]['deltaG'])
        y=np.array(df[df.paired == is_paired]['mut_rates'])

        tra[prefix] = go.Scatter(
            x=x,
            y=y,
            text = df[hover_attr],
            mode='markers',
            name=prefix,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)]),
            line=dict(color='green' if is_paired else 'red'))
            
        
        for m in models:
            if len(y) > 0:
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
        plot(fig, filename = savefile, auto_open=auto_open)
    return OutputPlot(df, fig)
    
def variable_exp_across_samples(df:pd.DataFrame, experimental_variable:str, models:List[str]=[], savefile=None, auto_open=False)->OutputPlot:

    data = pd.DataFrame()
    for _, row in df.iterrows():
        data = pd.concat([data, pd.DataFrame({'sample':row['sample'],experimental_variable:getattr(row,experimental_variable), 'index':list(row.filtered_index), 'base':list(row.sequence), 'mut_rates':list(row.mut_rates), 'paired':[s !='.' for s in row.structure_selected]}, index=list(range(len(row.filtered_index))))])
    data = data.reset_index().rename(columns={'level_0':'index_subsequence'})
    data = data.sort_values(by='index')
    data['Marker'] = data['paired'].apply(lambda x: {True: 0, False:1}[x])
    hover_attr = ['base','mut_rates','sample',experimental_variable]

    tra = {}
    for idx, row in data.groupby('index_subsequence'):
        name = f"({row['base'].iloc[0]},{row['index'].iloc[0]})"
        tra[row['index'].iloc[0]] = go.Scatter(
            x= row[experimental_variable], 
            y= row['mut_rates'], 
            text = data[hover_attr],
            mode='lines+markers',
            marker = dict(symbol = list(map(util.Setshape, data['Marker']))),
            name= name,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)])
        )         
        for m in models:
            x= row[experimental_variable]
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

    #tra = {k:tra[k] for k in sorted(tra.keys())}

    fig = dict(data = list(tra.values()), layout = layout)
    iplot(fig)
    if savefile != None:
        plot(fig, filename = savefile, auto_open=auto_open)
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

