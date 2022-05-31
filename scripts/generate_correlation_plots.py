import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from os.path import exists, dirname
import os, sys
from scipy.stats import linregress
from nap import *
try:
    sys.path.append(dirname('libs/dreem/dreem')) 
except:
    "If dreem isn't installed on your computer, the code won't run"
    
from matplotlib.offsetbox import AnchoredText

# Set your username for the database (at the moment, keep Yves)
username = 'Yves'

# Select your study
study = 'all_tubes' 

## Set your base coverage high-pass filter value
min_bases_cov = 1000 

# Set the resolution for the plots
mpl.rcParams['figure.dpi'] = 200 # the highest the resolution, the slowest the plotting

# Depending on the study you select, you'll get a series of tubes. You can also create new studies using this dictionary.
# Here's an example.
tubes_per_study = {   
    'replicates':           ['C5', 'A4', 'F4', 'A6', 'A7'],
    'salt':                 ['A6', 'B6', 'C6', 'D6', 'E6'], 
    'temperature':          ['D7', 'E7', 'F7', 'G7', 'H7', 'A8', 'B8', 'C8'], 
    'magnesium':            ['F6', 'G6', 'H6', 'A7', 'B7', 'C7'],
    '60 mM DMS kinestics':  ['D8', 'E8', 'F8', 'G8', 'H8', 'A9'],
    'all_tubes': [ele for ele in [f"{a}{b}" for a in string.ascii_uppercase[0:8] for b in range(1,11)] \
                 if ele not in ['C3','C10','D10','E10','F10','G10','H10', 'E4']]
                 + ['C5_realignment_v3']
    }

tubes = tubes_per_study[study]


# Pull the firebase
df_rough = firebase.load(tubes=tubes, username=username)

# Clean and reformat the dataset
df, df_full = data_wrangler.clean_dataset(df_rough=df_rough,
                                             tubes=tubes, 
                                             min_bases_cov=min_bases_cov)
                                             
                                   


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
            get_roi_info(df, tubes[1], construct)
            for roi in roi_structure_comparison_color:
                try:
                    x, y = np.array(get_roi_info(df, tubes[1], construct)['mut_rate'].xs((is_paired,roi), level=('paired','roi_structure_comparison')), dtype=float),\
                            np.array(get_roi_info(df, tubes[0], construct)['mut_rate'].xs((is_paired,roi),level=('paired','roi_structure_comparison')), dtype=float)
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
    fig, axs = plt.subplots(len(tubes)+1,len(tubes), figsize= (25,25), sharex=True, sharey=True)
    for x in range(1,len(tubes)+1):
        for y in range(0,len(tubes)):
            df_global_corr = pd.concat((df_global_corr,correlation_2_tubes(df, (tubes[x-1], tubes[y]), constructs, axs[x][y])),
                                        axis = 0,
                                        join="outer",
                                        ignore_index=True)
    axs[0,len(tubes)-2].plot(0,0,'b.',0,0,'r.',0,0,'bx',0,0,'rx',0,0,'g-')            
    axs[0,len(tubes)-2].legend(['Paired in full sequence RNAstructure, paired in ROI RNAstructure',
                'Paired in full sequence RNAstructure, not paired in ROI RNAstructure',
                'Not paired in full sequence RNAstructure, not paired in ROI RNAstructure',
                'Not paired in full sequence RNAstructure, paired in ROI RNAstructure',
                'Fit'])
    return df_global_corr
    
    
show_plots = False
for study in tubes_per_study:
    tubes = tubes_per_study[study]
    for constructs in df.construct.unique():
        constructs_name = constructs
        df_global_corr = correlation_n_tubes(df, tubes, constructs)
        plt.title(f"Correlation between tubes for study: {study}, constructs: {constructs_name}")
        plot.save_fig(path=f"data/figs/date/correlation/{study}/{constructs_name}", 
                        title=f"correlation_fit_{study}_{constructs_name}")
        plt.close(not show_plots)

        for plt_type in ['r_value', 'slope']:
            pivot = df_global_corr.pivot("tube_0","tube_1", plt_type).astype(float)
            f, ax = plt.subplots(figsize=(28, 10))
            sns.heatmap(pivot, annot=False, linewidths=0, ax=ax)#, norm=LinNorm())
            plt.title(f"{plt_type} of the correlation between tubes for study: {study}, constructs: {constructs_name}")
            plot.save_fig(path=f"data/figs/date/correlation/{study}/{constructs_name}", 
                            title=f"correlation_{plt_type}_{study}_{constructs_name}")
            plt.close(not show_plots)
              
