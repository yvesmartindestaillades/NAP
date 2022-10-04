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


def mut_histogram(df, samp:str, construct:str, gene:str=None, cluster:int=None, structure:str=None, show_ci:bool=True)->OutputPlot:

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

    plotly.offline.plot(
            mut_fig, filename= "pop_avg.html", auto_open=True,
    )
    return OutputPlot(mh, mut_fig)

