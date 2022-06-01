from matplotlib import colors
import pandas as pd
import pickle
import json
import firebase_admin
from firebase_admin import credentials
from firebase_admin import db
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import string
from os.path import exists
import os
import datetime
import seaborn as sns
from os.path import exists, dirname
import os, sys

try:
    sys.path.append(dirname('/libs/dreem/dreem')) 
except:
    "If dreem isn't installed on your computer, the code won't run"


from libs import dreem
from scipy.stats import linregress
from matplotlib.offsetbox import AnchoredText
sys.path.append(os.path.abspath(""))
from nap import utils


def tube_coverage_distribution(df):
    plt.figure()
    plt.plot(np.array(df['tubes_covered'].sort_values(ascending=False)))
    plt.xlabel('Constructs (sorted)')
    plt.ylabel('Amount of tubes covered')
    plt.grid()

def valid_construct_per_tube(df, min_bases_cov):
    df.groupby('tube').count().reset_index().plot(kind='bar',x='tube', y='construct', figsize=figsize )
    plt.ylabel(f"Number of construct above {min_bases_cov} reads")
    plt.grid()

def base_coverage_for_all_constructs(df, min_bases_cov):
    plt.figure()
    plt.plot(np.array(df['cov_bases_roi'].sort_values(ascending=False).reset_index())[:,1])
    plt.plot(np.arange(0,int(df.construct.count()),1), [min_bases_cov]*int(df.construct.count()))
    plt.legend(['Dataframe', '1000 reads line'])
    plt.xlabel('Constructs (sorted)')
    plt.ylabel('# of reads of the worst covered base in the ROI for a given structure in a given tube')

def random_9_base_coverage(df, min_bases_cov):
    random_selection = np.random.randint(len(df), size=(9))
    fig = plt.figure()
    for i in range(9):
        axes1 = plt.subplot(int('33'+str(i+1)))
        plt.plot(np.array(df['cov_bases'].iloc[random_selection[i]]))
        start, end = df['roi_start_index'].iloc[random_selection[i]], df['roi_end_index'].iloc[random_selection[i]]
        plt.plot(np.arange(start, end, 1), np.array(df['cov_bases'].iloc[random_selection[i]])[start:end])
        plt.plot(np.arange(0, len(df['cov_bases'].iloc[random_selection[i]])), len(df['cov_bases'].iloc[random_selection[i]])*[min_bases_cov])
        plt.xlabel("Bases")
        plt.ylabel("Coverage")
        plt.title(f"Construct {df['construct'].iloc[random_selection[i]]}, tube {df['tube'].iloc[random_selection[i]]} ")
        plt.grid()
        plt.legend(["Base coverage (all)", 'Base coverage (ROI)', 'min_bases_cov'])
        axes2 = axes1.twinx()   
        axes2.set_ylabel('Coverage [%]')
        axes2.set_ylim((0,100*max(df['cov_bases'].iloc[random_selection[i]]))/df['num_reads'].iloc[random_selection[i]])
    fig.tight_layout()

def base_coverage(df, tube, construct, min_bases_cov=None):
    ax1 = plt.subplot()
    serie = df.set_index(['tube','construct']).loc[tube, construct]
    plt.plot(np.array(serie['cov_bases']))
    start, end = serie['roi_start_index'], serie['roi_end_index']
    plt.plot(np.arange(start, end, 1), np.array(serie['cov_bases'])[start:end])
    if min_bases_cov != None:
        plt.plot(np.arange(0, len(serie['cov_bases'])), len(serie['cov_bases'])*[min_bases_cov])
    plt.xlabel("Bases")
    plt.ylabel("Coverage")
    plt.title(f"Construct {construct}, tube {tube} ")
    plt.grid()
    plt.legend(["Base coverage (all)", 'Base coverage (ROI)', 'min_bases_cov'])
    ax2 = ax1.twinx()   
    ax2.set_ylabel('Coverage [%]')
    ax2.set_ylim(0,100*max(serie['cov_bases'])/serie['num_reads'])
    plt.tight_layout()


def heatmap(df, column):
    base_cov_plot = df.pivot("tube","construct", column).astype(float)
    f, ax = plt.subplots()
    sns.heatmap(base_cov_plot, annot=False, linewidths=0, ax=ax, norm=LogNorm())


def fit_deltaG(df, tube): #TODO
            # Fit
    fit = lambda a, b, c, dG: a/(1+b*np.exp(-dG/(R*T))) + c

        ## Max mut freq
    #a = np.array(df_use['mod_bases_A'].loc[tube][-1]+df_use['mod_bases_C'].loc[tube][-1]).mean()/np.array(df_use['info_bases'].loc[tube][-1]).mean()
        ## Baseline mut freq
    #   c = np.array(df_use['mod_bases_A'].loc[tube][0]+df_use['mod_bases_C'].loc[tube][0]).mean()/np.array(df_use['info_bases'].loc[tube][0]).mean()
        ## #TODO b
    # b = 1
    #  print(f"For tube {tube}, max mut freq a is {a}, b is {b}, baseline mut freq c is {c}")
    # plt.plot(t, fit(a, b, c, t))



def mutation_rate(df, tube, construct, plot_type, index):
    
    df_use = df.set_index(['tube','construct'])
    
    if not plot_type in ['sequence','partition']:
        raise f"{plot_type} must be 'sequence' or 'partition', please check this argument"

    if plot_type == 'sequence':  # Plot the mutation rate for each base along the sequence

        mut_per_base = pd.DataFrame({'mut_rate': pd.Series(np.array(df_use[f"mut_bases"].loc[tube, construct][1:])/np.array(df_use[f"info_bases"].loc[tube, construct][1:]), dtype=object)
                                    ,'base':list(df_use['full_sequence'].loc[tube, construct])})\
                                    .reset_index()\
                                    .set_index(['base', 'index'])

        df_hist = pd.DataFrame()
        df_hist.index = mut_per_base.reset_index()['index']

        for base in ['A','C','G','T']:
            df_hist[base] = pd.Series(dtype=float)
            df_hist[base] = mut_per_base.loc[base]

        if index == 'base':
            df_hist.index = mut_per_base.reset_index()['base']

        ax = df_hist.plot.bar(stacked=True, figsize=(35,7), color=['r','b','y','g'])
        plt.title(f"tube {tube}, construct {construct}")

    if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
        df_hist = pd.DataFrame()
        for base in ['A','C','G','T']:
            df_hist[f"mod_bases_{base}"]  = np.array(df_use[f"mod_bases_{base}"].loc[tube, construct][1:])/df_use['info_bases'].loc[tube, construct][1:]

        if index == 'base':
            df_hist.index = list(df_use['full_sequence'].loc[tube,construct])

        ax = df_hist.plot.bar(stacked=True, figsize=(35,7), color=['r','b','y','g'])



def deltaG(df, tube):
    df_use = df.set_index(['tube','construct'])

    fig = define_figure(title=tube,
                            xlabel='deltaG [kcal]',
                            ylabel='Mutation ratio',
                            figsize=(20,5))

    stack_for_plot = {'0':{'x':[],'y':[]},'1':{'x':[],'y':[]}}

    for construct in df.construct.unique():
        roi_part = utils.get_roi_info(df=df, tube=tube, construct=construct)
        for base in ['A','C']:
            for roi_struct_comp in ['0','1']:
                try:    
                    this_base_mut =  roi_part.xs((base,True,roi_struct_comp), level=('base','paired','roi_structure_comparison'))
                    stack_for_plot[roi_struct_comp]['x'].extend(this_base_mut['roi_deltaG'].to_list())
                    stack_for_plot[roi_struct_comp]['y'].extend(this_base_mut['mut_rate'].to_list())
                except:
                    continue
    plt.plot(stack_for_plot['0']['x'],stack_for_plot['0']['y'],'b.')
    plt.plot(stack_for_plot['1']['x'],stack_for_plot['1']['y'],'r.')
    fit_deltaG(df_use, tube)
    plt.legend(['A and C bases of the ROI, predicted paired by RNAstructure for both the ROI sequence and the full sequence',\
                'A and C bases of the ROI part, predicted paired by RNAstructure for the full sequence but not for the ROI sequence'])
    plt.ylim([0,0.15])
    fig.tight_layout()


def correlation_2_tubes(df, tubes, constructs, axs=None):

    if type(constructs) != list:
        constructs = [constructs]

    if axs is None:
        fig, axs = plt.subplots(1,1)

    paired = {True: '.',False:'x'}
    roi_structure_comparison_color = {'0':'b','1':'r'}
    x_all, y_all = [], []
    for construct in constructs:
        for is_paired in paired: 
            utils.get_roi_info(df, tubes[1], construct)
            for roi in roi_structure_comparison_color:
                try:
                    x, y = np.array(utils.get_roi_info(df, tubes[1], construct)['mut_rate'].xs((is_paired,roi), level=('paired','roi_structure_comparison')), dtype=float),\
                            np.array(utils.get_roi_info(df, tubes[0], construct)['mut_rate'].xs((is_paired,roi),level=('paired','roi_structure_comparison')), dtype=float)
                    axs.plot(x,y,f"{roi_structure_comparison_color[roi]}{paired[is_paired]}")
                    axs.tick_params(axis='x', labelrotation = 45)
                    x_all.extend(x), y_all.extend(y)
                except:
                    axs.plot()
                    continue
    result = linregress(x_all,y_all)
    p =  np.poly1d((result.slope,result.intercept))
    t = np.linspace(min(x_all),max(x_all))
    axs.plot(t,p(t),'g-')
    axs.grid()
    axs.set(xlabel=f"Mutation rate of tube {tubes[1]}", ylabel=f"Mutation rate of tube {tubes[0]}")
    anchored_text = AnchoredText(f"R = {round(result.rvalue,3)}, slope = {round(result.slope,3)}", loc=2)
    axs.add_artist(anchored_text)
    df_global_corr = pd.DataFrame({'tube_0':tubes[0], 'tube_1':tubes[1], 'r_value':result.rvalue, 'slope':result.slope}, index=[0])
    return df_global_corr

def correlation_n_tubes(df, tubes, constructs):  
    df_global_corr = pd.DataFrame(columns=['tube_0', 'tube_1', 'r_value', 'slope'])
    plt.figure(figsize=(30*len(tubes), 30*len(tubes)))
    fig, axs = plt.subplots(len(tubes)+1,len(tubes), figsize= (30,30), sharex=True, sharey=True)
    for x in range(1,len(tubes)+1):
        for y in range(0,len(tubes)):
            df_global_corr = pd.concat((df_global_corr, correlation_2_tubes(df, (tubes[x-1], tubes[y]), constructs, axs[x][y])),
                                        axis = 0,
                                        join="outer",
                                        ignore_index=True)
    axs[0,len(tubes)-1].plot(0,0,'b.',0,0,'r.',0,0,'bx',0,0,'rx',0,0,'g-')            
    axs[0,len(tubes)-1].legend(['Paired in full sequence RNAstructure, paired in ROI RNAstructure',
                'Paired in full sequence RNAstructure, not paired in ROI RNAstructure',
                'Not paired in full sequence RNAstructure, not paired in ROI RNAstructure',
                'Not paired in full sequence RNAstructure, paired in ROI RNAstructure',
                'Fit'])
    return df_global_corr

def save_fig(path,title):
    full_path = utils.make_path(path)
    plt.savefig(f"{full_path}/{title}")

def define_figure(title, xlabel, ylabel, figsize):
    fig = plt.figure(figsize=figsize)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return fig
