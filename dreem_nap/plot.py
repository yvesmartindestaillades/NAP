import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from os.path import exists, dirname
import os, sys
from scipy.stats import linregress
from matplotlib.offsetbox import AnchoredText

from dreem_nap.study import Study
sys.path.append(os.path.abspath(""))
from dreem_nap import data_manip, data_wrangler, utils
from typing import Tuple, List



def sample_coverage_distribution(df:pd.DataFrame)->None:
    """Plot each construct vs its number of samples covered, i.e the number of samples containing this construct.
    
    Args:
        df: dataframe of interest.
    """

    plt.figure()
    plt.plot(np.array(df.groupby['construct']['samples_covered'].mean().sort_values(ascending=False)))
    plt.xlabel('Constructs (sorted)')
    plt.ylabel('Amount of samples covered')
    plt.grid()


def valid_construct_per_sample(df:pd.DataFrame)->None:
    """Plot how many valid constructs each sample has.

    Args:
        df: dataframe of interest.
    
    Raises:
        If the min_bases_cov value isn't the same for every samples.
    """
    
    df.groupby('samp').count().reset_index().plot(kind='bar',x='samp', y='construct')
    assert len(df.min_bases.cov.unique()) == 1, "min_bases_cov isn't unique"
    plt.ylabel(f"Number of construct above {int(df.min_bases_cov.unique())} reads")
    plt.grid()


def base_coverage_ROI_for_all_constructs(df:pd.DataFrame)->None:
    """Plot the base-coverage of the worst-covered base of the Region of Interest, for each construct. 
    
    Args:
        df: dataframe of interest.
    """
    plt.figure()
    plt.plot(np.array(df['cov_bases_roi'].sort_values(ascending=False).reset_index())[:,1])
    plt.plot(np.arange(0,int(df.construct.count()),1), [int(df.min_bases_cov.unique())]*int(df.construct.count()))
    plt.legend(['Dataframe', '1000 reads line'])
    plt.xlabel('Constructs (sorted)')
    plt.ylabel('# of reads of the worst covered base in the ROI for each (sample, construct)')


def random_9_base_coverage(df:pd.DataFrame)->None:
    """Plot the base coverage of 9 randomly picked (sample, construct).

    Args:
        df: dataframe of interest.
    """

    random_selection = np.random.randint(len(df), size=(9))
    fig = plt.figure(figsize=(25,10))
    for i in range(9):
        axes1 = plt.subplot(int('33'+str(i+1)))
        plt.plot(np.array(df['cov_bases'].iloc[random_selection[i]]))
        start, end = df['roi_start_index'].iloc[random_selection[i]], df['roi_end_index'].iloc[random_selection[i]]
        plt.plot(np.arange(start, end, 1), np.array(df['cov_bases'].iloc[random_selection[i]])[start:end])
        plt.plot(np.arange(0, len(df['cov_bases'].iloc[random_selection[i]])), len(df['cov_bases'].iloc[random_selection[i]])*[int(df.min_bases_cov.unique())])
        plt.xlabel("Bases")
        plt.ylabel("Coverage")
        plt.title(f"Construct {df['construct'].iloc[random_selection[i]]}, sample {df['samp'].iloc[random_selection[i]]} ")
        plt.grid()
        plt.legend(["Base coverage (all)", 'Base coverage (ROI)', 'min_bases_cov'])
        axes2 = axes1.twinx()   
        axes2.set_ylabel('Coverage [%]')
        axes2.set_ylim((0,100*max(df['cov_bases'].iloc[random_selection[i]]))/df['num_reads'].iloc[random_selection[i]])
    fig.tight_layout()


def base_coverage(df:pd.DataFrame, samp:str, construct:int)->None:
    """Plot the base coverage of a specific (sample, construct).
    
    Args:
        df: dataframe of interest.
        samp: sample of interest.
        construct: construct of interest.   
    """

    ax1 = plt.subplot()
    serie = df.set_index(['samp','construct']).loc[samp, construct]
    plt.plot(np.array(serie['cov_bases']))
    start, end = serie['roi_start_index'], serie['roi_end_index']
    plt.plot(np.arange(start, end, 1), np.array(serie['cov_bases'])[start:end])
    plt.plot(np.arange(0, len(serie['cov_bases'])), len(serie['cov_bases'])*[int(df.min_bases_cov.unique())])
    plt.xlabel("Bases")
    plt.ylabel("Coverage")
    plt.title(f"Construct {construct}, sample {samp} ")
    plt.grid()
    plt.legend(["Base coverage (all)", 'Base coverage (ROI)', 'min_bases_cov'])
    ax2 = ax1.twinx()   
    ax2.set_ylabel('Coverage [%]')
    ax2.set_ylim(0,100*max(serie['cov_bases'])/serie['num_reads'])
    plt.tight_layout()


def heatmap(df:pd.DataFrame, column:str)->None:
    """Plot the heatmap of a specific column of your dataframe, w.r.t the samples and constructs.

    Args:
        df: a Pandas dataframe.
        column: a specific column of your dataframe.    
    """

    base_cov_plot = df.pivot("samp","construct", column).astype(float)
    f, ax = plt.subplots()
    sns.heatmap(base_cov_plot, annot=False, linewidths=0, ax=ax, norm=LogNorm())


def mut_rate_vs_base_non_pairing_prob(df:pd.DataFrame, samp:str, construct:int, plot_type=['mut','prob','corr'], bases_type=['A','C','G','T'], roi_range = list(range(0,170))):
    """Plot a mutation rate histogram, a base non-pairing probability histogram, and a scatter plot fitting the mutation rate vs the base non-pairing. 
       The base non-pairing values are 1 - base-pairing values.

    Args:
        df (pd.DataFrame): your dataset
        samp (str): your sample of interest
        construct (int): your construct of interest.
        plot_type (list[str]): which plots that you want. Default is ['mut','prob','corr'] (all plots). \
        'mut' is for the mutation histogram. 'prob' is for the non-pairing probability histogram. 'corr' is for the scatter plot correlating \
        the mutation rate and the non-pairing probability. 'merge' gives all plots together.
        bases_type (list[str]): which bases that you want to take into account. Default is ['A','C','G','T']
        roi_range (list[int]): 0-index of the bases that you want to take into account. Default is [0, 1, ..., 168, 169].
    """

    COLOR_PER_BASE = {'A':'r','C':'b','G':'y','T':'g'}

    def split_bases(df, column):
        """Split a sample-construct dataframe into 4 columns, one per base. The columns contain the value of the column 'column' entered as an input.
        """

        df = pd.DataFrame(df[column])
        df_hist = pd.DataFrame()
        df_hist.index = df.reset_index()['index']

        for base in bases_type:
            df_hist[base] = pd.Series(dtype=float)
            df_hist[base] = df.loc[base]

        return df_hist

    def plot_histogram(df, title, ax = None):
        """Plot histogram using a base-wise splitted sample-construct dataframe.
        """

        df.plot.bar(stacked=True, figsize=(35,7), color= [COLOR_PER_BASE[base] for base in bases_type], rot=90, ax=ax)
        plt.title(f"{title}, sample {samp}, construct {construct}")
    
    def plot_scatter(df):
        """Plot scatter plot using a base-wise splitted sample-construct dataframe.
        """
        
        from sklearn.linear_model import LinearRegression

        x, y = df['base non-pairing probability'], df['mut_rate']
        x, y = np.array(x).reshape(-1,1), np.array(y).reshape(-1,1)

        for base, color in zip(bases_type, [COLOR_PER_BASE[base] for base in bases_type]):
            px, py = split_bases(df, 'base non-pairing probability'), split_bases(df, 'mut_rate')
            plt.plot(px[base], py[base],color+'.')
        plt.legend(['A','C','G','T'])
        plt.xlabel('Base non-pairing probability')
        plt.ylabel('Mutation rate')

        model = LinearRegression().fit(x,y)
        plt.plot(model.predict(np.array([0,1]).reshape(-1,1)))
        plt.title(f"Base non-pairing prob vs mut rate \n \
    sample {samp}, construct {construct} \n " + \
    'R2= {:.5}, y = {:.5} x + {:.5}'.format(model.score(x, y), float(model.coef_[0]),model.intercept_[0])) 


    # Data wrangling
    df_sc = data_manip.get_roi_info(df, samp, construct, bases_type=bases_type, roi_range=roi_range)\
                                    .reset_index()\
                                    .drop(columns=['paired','roi_structure_comparison','roi_deltaG'])\
                                    .set_index(['base','index'])
    df_sc['base non-pairing probability'] = 1 - df_sc['base_pairing_prob']

    # Plots
    if 'prob' in plot_type:
        plot_histogram(split_bases(df_sc, 'base non-pairing probability'), 'Base non-pairing probability')
    if 'mut' in plot_type:
        plot_histogram(split_bases(df_sc, 'mut_rate'), 'Mutation rate')
    if 'corr' in plot_type:
        plot_scatter(df_sc)
    
    if 'merge' in plot_type:
        plt.figure(figsize=(20,10))
        plot_histogram(split_bases(df_sc, 'base non-pairing probability'), 'Base non-pairing probability', plt.subplot(2,1,1))
        plot_histogram(split_bases(df_sc, 'mut_rate'), 'Mutation rate', plt.subplot(2,1,2))
        plt.tight_layout()


def mut_histogram(plot_type:str, df:pd.DataFrame, samp:str, construct:int)->None:
    """Plot the mutation rate of a specific (sample, construct).

    Args:
        plot_type: 'sequence' or 'partition'. 
            - 'index' uses bases numbers as index and the original construct bases as colors.
            - 'partition' uses original sequence bases as index and the partition of mutated bases as colors.
        df: dataframe of interest.
        samp]: sample of interest.
        construct: construct of interest.
    """

    df_use = df.set_index(['samp','construct'])
    
    if not plot_type in ['index','partition']:
        raise f"{plot_type} must be 'index' or 'partition', please check this argument"

    if plot_type == 'index':  # Plot the mutation rate for each base along the sequence

        mut_per_base = pd.DataFrame({'mut_rate': pd.Series(np.array(df_use[f"mut_bases"].loc[samp, construct][1:])/np.array(df_use[f"info_bases"].loc[samp, construct][1:]), dtype=object)
                                    ,'base':list(df_use['full_sequence'].loc[samp, construct])})\
                                    .reset_index()\
                                    .set_index(['base', 'index'])

        df_hist = pd.DataFrame()
        df_hist.index = mut_per_base.reset_index()['index']

        for base in ['A','C','G','T']:
            df_hist[base] = pd.Series(dtype=float)
            df_hist[base] = mut_per_base.loc[base]

        df_hist.index = mut_per_base.reset_index()['base']

        ax = df_hist.plot.bar(stacked=True, figsize=(35,7), color=['r','b','y','g'])
        plt.title(f"sample {samp}, construct {construct}")

    if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
        df_hist = pd.DataFrame()
        for base in ['A','C','G','T']:
            df_hist[f"mod_bases_{base}"]  = np.array(df_use[f"mod_bases_{base}"].loc[samp, construct][1:])/df_use['info_bases'].loc[samp, construct][1:]

        df_hist.index = list(df_use['full_sequence'].loc[samp,construct])

        ax = df_hist.plot.bar(stacked=True, figsize=(35,7), color=['r','b','y','g'])


def deltaG(df:pd.DataFrame, samp:str, bases_type=['A','C'])->None:
    """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.
    
    Args:
        df: dataframe of interest.
        sample: sample of interest.
    """

    fig = utils.define_figure(title=samp,
                            xlabel='deltaG [kcal]',
                            ylabel='Mutation ratio',
                            figsize=(20,5))

    stack_for_plot = {'0':{'x':[],'y':[]},'1':{'x':[],'y':[]}}

    for construct in df[df.samp==samp].construct.unique():
        roi_part = data_manip.get_roi_info(df=df, samp=samp, construct=construct, bases_type=bases_type)
        for base in bases_type:
            for roi_struct_comp in ['0','1']:
                try:    
                    this_base_mut =  roi_part.xs((base,True,roi_struct_comp), level=('base','paired','roi_structure_comparison'))
                    stack_for_plot[roi_struct_comp]['x'].extend(this_base_mut['roi_deltaG'].to_list())
                    stack_for_plot[roi_struct_comp]['y'].extend(this_base_mut['mut_rate'].to_list())
                except:
                    continue
    plt.plot(stack_for_plot['0']['x'],stack_for_plot['0']['y'],'b.')
    plt.plot(stack_for_plot['1']['x'],stack_for_plot['1']['y'],'r.')
    plt.legend(['A and C bases of the ROI, predicted paired by RNAstructure for both the ROI sequence and the full sequence',\
                'A and C bases of the ROI part, predicted paired by RNAstructure for the full sequence but not for the ROI sequence'])
    plt.ylim([0,0.15])
    fig.tight_layout()

    
def deltaG_basewise(df:pd.DataFrame, samp:str, roi_range, sublib = None, full_seq_pairing=None, bases_type = ['A','C'], style='.', figsize=(20,5))->None:
    """Plot the mutation rate of each paired-predicted base of the ROI for each construct of a sample, w.r.t the deltaG estimation.
    
    Args:
        df: dataframe of interest.
        sample: sample of interest.
        roi_range: 0-index index of the bases to plots
        sublib: a str contained in your sub-library(ies) of interest.
        full_seq_pairing: pairing state of plotted bases, w.r.t the full-seq prediction. default is None. None: both paired and unpaired bases. 'paired': only paired bases. 'unpaired': only unpaired bases.
    """

    # Define utils
    PAIRED, UNPAIRED = '1','0' 
    to_str = {PAIRED: 'paired', UNPAIRED:'unpaired'}  
    to_macro = {v: k for k, v in to_str.items()} 
    xor = lambda a,b: bool(((not a) and b) or ((not b) and a))

    # filter by sub-library
    if sublib != None:
        df = df.loc[df['sub-library'].apply(lambda x: sublib in x)]

    # Select which pairing
    sequence_pairings = [PAIRED, UNPAIRED]
    if full_seq_pairing != None:
        sequence_pairings = [to_macro[full_seq_pairing]]                   

    # Utils for plotting
    legend = []
    fig = utils.define_figure(title=samp,
                            xlabel='deltaG [kcal]',
                            ylabel='Mutation ratio',
                            figsize=figsize)

    stack_for_plot = {PAIRED: {PAIRED:{'x':[],'y':[]},UNPAIRED:{'x':[],'y':[]}}, 
                      UNPAIRED: {PAIRED:{'x':[],'y':[]},UNPAIRED:{'x':[],'y':[]}}}

    # For each construct of this sample get your data
    for construct in df[df.samp == samp].construct.unique():
        roi_part = data_manip.get_roi_info(df=df, samp=samp, construct=construct, roi_range=roi_range, bases_type=bases_type)
        for sequence_pairing in sequence_pairings:  # Sequence-wise prediction
            for roi_struct_comp in [UNPAIRED, PAIRED]: 
                roi_pairing = str(int( xor(bool(int(roi_struct_comp)), bool(int(sequence_pairing))))) # ROI-wise prediction
                try:    
                    this_base_mut = roi_part.xs((int(sequence_pairing),roi_struct_comp), level=('paired','roi_structure_comparison')).reset_index()
                    this_base_mut = this_base_mut[this_base_mut['index'].isin(roi_range)]
                    stack_for_plot[sequence_pairing][roi_pairing]['x'].extend(this_base_mut['roi_deltaG'].to_list())
                    stack_for_plot[sequence_pairing][roi_pairing]['y'].extend(this_base_mut['mut_rate'].to_list())
                except:
                    continue

    # Plot the data
    for sequence_pairing in sequence_pairings:  # Sequence-wise prediction
        for roi_struct_comp in [UNPAIRED, PAIRED]: 
            roi_pairing = str(int( xor(bool(int(roi_struct_comp)), bool(int(sequence_pairing))))) # ROI-wise prediction
            plt.plot(stack_for_plot[sequence_pairing][roi_pairing]['x'],stack_for_plot[sequence_pairing][roi_pairing]['y'], style)
            legend += [f"A and C selected bases, sequence-wise pairing prediction: {to_str[sequence_pairing]}, ROI-wise pairing prediction: sequence-wise pairing prediction: {to_str[roi_pairing]}"]
    plt.legend(legend)
    plt.ylim([0,0.15])
    fig.tight_layout()




def _correlation_2_samples(df:pd.DataFrame, samples:Tuple[str,str], constructs:List[int], axs:plt.axes=None)->pd.DataFrame:
    """Plot the mutation rate of each paired-predicted base of the ROI for a sample vs the same base in another sample, and this for specific constructs.

    Args:
        df: dataframe of interest.
        samples: samples of interest.
        constructs: constructs of interest.
        axs: a plt.axes object to use - for example for subplotting.
    
    Returns:
        A single-lined Pandas dataframe with values for r_value and slope.
    """

    if axs is None:
        fig, axs = plt.subplots(1,1)

    paired = {True: '.',False:'x'}
    roi_structure_comparison_color = {'0':'b','1':'r'}
    x_all, y_all = [], []
    for construct in constructs:
        for is_paired in paired: 
            data_manip.get_roi_info(df, samples[1], construct)
            for roi in roi_structure_comparison_color:
                try:
                    x, y = np.array(data_manip.get_roi_info(df, samples[1], construct)['mut_rate'].xs((is_paired,roi), level=('paired','roi_structure_comparison')), dtype=float),\
                            np.array(data_manip.get_roi_info(df, samples[0], construct)['mut_rate'].xs((is_paired,roi),level=('paired','roi_structure_comparison')), dtype=float)
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
    axs.set(xlabel=f"Mutation rate of sample {samples[1]}", ylabel=f"Mutation rate of sample {samples[0]}")
    anchored_text = AnchoredText(f"R = {round(result.rvalue,3)}, slope = {round(result.slope,3)}", loc=2)
    axs.add_artist(anchored_text)
    df_corr = pd.DataFrame({'sample_0':samples[0], 'sample_1':samples[1], 'r_value':result.rvalue, 'slope':result.slope}, index=[0])
    return df_corr


def correlation_n_samples(df:pd.DataFrame, study:Study, constructs:List[int])->pd.DataFrame:  
    """Plot correlation_2_samples() for each possible pair in the samples list, and this for specific constructs, in a single plot.

    Args:
        df: dataframe of interest.
        study: class containing a list of samples that you want to use.
        constructs: constructs of interest.

    Returns:
        A Pandas dataframe with values for r_value and slope for each possible pair in the samples list.
    """
    samples = study.samples

    df_global_corr = pd.DataFrame(columns=['sample_0', 'sample_1', 'r_value', 'slope'])
    plt.figure(figsize=(30*len(samples), 30*len(samples)))
    fig, axs = plt.subplots(len(samples)+1,len(samples), figsize= (30,30), sharex=True, sharey=True)
    for x in range(1,len(samples)+1):
        for y in range(0,len(samples)):
            df_global_corr = pd.concat((df_global_corr, _correlation_2_samples(df, (samples[x-1], samples[y]), constructs, axs[x][y])),
                                        axis = 0,
                                        join="outer",
                                        ignore_index=True)
    axs[0,len(samples)-1].plot(0,0,'b.',0,0,'r.',0,0,'bx',0,0,'rx',0,0,'g-')            
    axs[0,len(samples)-1].legend(['Paired in full sequence RNAstructure, paired in ROI RNAstructure',
                'Paired in full sequence RNAstructure, not paired in ROI RNAstructure',
                'Not paired in full sequence RNAstructure, not paired in ROI RNAstructure',
                'Not paired in full sequence RNAstructure, paired in ROI RNAstructure',
                'Fit'])
    return df_global_corr


def study_sample(df:pd.DataFrame, study:Study, bases_type:list[str]=['A','C'], sub_lib:str=None, structure:str='full', scale_x:str = 'lin', roi_range=None, overlay=0, figsize:Tuple[float,float]=None):
    """Plot the mean of the mutation rate of the ROI bases, for each sample of the study.

    Args:
        df (pd.DataFrame): dataframe of interest.
        study (Study): class containing relevant information about the series of sample that you want to use.
        bases_type (list[str], optional): list of the bases to filter-in. Defaults to ['A','C'].
        sub_lib (str, optional): select a specific set of constructs by grouping by the 'sub_library' column of the df. Defaults to None.
        structure (str, optional): what sequence the structure prediction is based on. Can be 'full' or 'roi'. Defaults to 'full'.
        roi_range: default is None. Array of base indexes (list[int]). ex: [80, 83, 91].
        scale_x (str, optional):  linear 'lin' or log 'log' for the x scale of plots. Defaults to 'lin'.
        overlay (int, optional): extend/shrink roi. Defaults to 0.
            'all': the roi is all bases
            int-type argument: the roi is the subsequence [start_roi_index-overlay, end_roi_index+overlay] 
            tuple[int]-type argument: the roi is the subsequence [start_roi_index+overlay[0], end_roi_index+overlay[1]].
        figsize (Tuple, optional): size of the plotted figure. Defaults to None.

    Raises:
        f: if sub-library isn't a valid element of the 'sub-library' column of the dataframe.
    """

    if sub_lib == None:
        constructs = df.construct.unique()
    else:
        if not sub_lib in df['sub-library'].unique():
            raise f"arg {sub_lib} is not a sub-library of the df"

        constructs = df[df['sub-library'] ==sub_lib].construct.unique()

    roi_range, roi_range_name = utils.roi_range_calc(overlay, roi_range, 
                                                roi_bounds=[df[df.construct==construct]['roi_start_index'].iloc[0], df[df.construct==construct]['roi_end_index'].iloc[0]],
                                                full_bounds=[df[df.construct==construct]['start'].iloc[0]-1, df[df.construct==construct]['end'].iloc[0]-1])

    samples = study.samples
    paired, unpaired = np.zeros(len(samples)), np.zeros(len(samples))
    for count, samp in enumerate(samples):
        for construct in constructs:
            df_roi = data_manip.get_roi_info(df, samp, construct, bases_type=bases_type,structure=structure, roi_range=roi_range)
            if True in df_roi.reset_index()['paired'].unique():
                paired[count] += df_roi.xs(True, level= 'paired')['mut_rate'].mean()/len( df.construct.unique())
            if False in df_roi.reset_index()['paired'].unique():
                unpaired[count] += df_roi.xs(False, level= 'paired')['mut_rate'].mean()/len( df.construct.unique())
    if figsize != None:
        fig = plt.figure(figsize=figsize)
    else:
        fig = plt.figure()
    
    scaled_plot = {'lin':plt.plot, 'log':plt.semilogx}[scale_x]
    scaled_plot(study.conditions, paired, 'r.-')
    scaled_plot(study.conditions, unpaired, 'b.-')
    plt.xlabel(f"{study.title}")
    plt.ylabel('Mean mutation rate for the ROI')
    plt.legend(['Paired-predicted bases','Unpaired-predicted bases'])

        
    fig.suptitle(f"Study: {study.name}, sub-library: {sub_lib}, bases types: {bases_type}, structure prediction: {structure}, roi_range: {roi_range_name}", fontsize=16)






def study_base(df:pd.DataFrame, study:Study, construct:int, bases_type = ['A','C'], structure = 'roi', scale_x = 'lin', roi_range:List[int]=None, overlay=0, figsize=(24,10))->None:
    """Generate line-plots of each base's mutation rate w.r.t a study's conditions, for a specific construct.

    Args:
        df (pd.DataFrame): dataframe of interest.
        study (Study): class containing relevant information about the series of sample that you want to use.
        construct (int): construct of interest.
        bases_type (list[str]): bases to display, sublist of ['A','C','G','T']
        structure: what sequence the structure prediction is based on. Can be 'full' or 'roi'.
        scale_x (str): linear 'lin' or log 'log'
        roi_range: default is None. Array of base indexes (list[int]). ex: [80, 83, 91].
        overlay (str or int or tuple[int]): extend/shrink roi
            'all': the roi is all bases
            int-type argument: the roi is the subsequence [start_roi_index-overlay, end_roi_index+overlay] 
            tuple[int]-type argument: the roi is the subsequence [start_roi_index+overlay[0], end_roi_index+overlay[1]].
        figsize (Tuple(int,int)): size of the plotted figure.
    
    Raise:
        "overlay and roi_range are non-compatible arguments."
    """

    assert not (overlay != 0 and roi_range != None), "overlay and roi_range are non-compatible arguments."

    df_paired, df_not_paired = pd.DataFrame(), pd.DataFrame()

    roi_range, roi_range_name = utils.roi_range_calc(overlay, roi_range, 
                                                roi_bounds=[df[df.construct==construct]['roi_start_index'].iloc[0], df[df.construct==construct]['roi_end_index'].iloc[0]],
                                                full_bounds=[df[df.construct==construct]['start'].iloc[0]-1, df[df.construct==construct]['end'].iloc[0]-1])

    for samp in study.samples:
        df_roi = data_manip.get_roi_info(df=df, samp=samp, construct=construct, bases_type=bases_type, structure=structure, roi_range=roi_range)
        if True in df_roi.reset_index()['paired'].unique():
            df_paired = pd.concat((df_paired, 
                              df_roi['mut_rate'].xs(True, level='paired').reset_index().set_index('index')
                              .drop(columns=['base','roi_structure_comparison']).transpose()))
        if False in df_roi.reset_index()['paired'].unique():
            df_not_paired = pd.concat((df_not_paired, 
                                   df_roi['mut_rate'].xs(False, level='paired').reset_index().set_index('index')
                                   .drop(columns=['base','roi_structure_comparison']).transpose()))


    intersection = lambda lst1, lst2: list(set(lst1) & set(lst2))

    if not df_paired.empty:
        df_paired = df_paired.set_index(pd.Series(study.conditions))
        paired_idx = intersection(df_paired.columns, roi_range)
        df_paired = df_paired[paired_idx].sort_index(axis=1)

    if not df_not_paired.empty:
        df_not_paired = df_not_paired.set_index(pd.Series(study.conditions))
        not_paired_idx = intersection(df_not_paired.columns, roi_range)    
        df_not_paired = df_not_paired[not_paired_idx].sort_index(axis=1)

    # Plot it
    fig = plt.figure()
    fig.suptitle(f"Construct: {construct}, study: {study.name}, bases types: {bases_type}, structure prediction: {structure}, roi_range: {roi_range_name}", fontsize=16)
    ax1, ax2 = plt.subplot(121), plt.subplot(122)

    if not df_paired.empty:
        df_paired.plot(figsize=figsize,
                        logx={'lin':False, 'log':True}[scale_x],
                        ax=ax1, 
                        title='Paired bases',  
                        xlabel=f"{study.title}",
                        ylabel="Mutation rate", 
                        style='.-')
    
    if not df_not_paired.empty:
        df_not_paired.plot(figsize=figsize,
                        logx={'lin':False, 'log':True}[scale_x],
                        ax=ax2, 
                        title='Unpaired bases', 
                        xlabel=f"{study.title}",
                        ylabel="Mutation rate",
                        style='.-')
    
