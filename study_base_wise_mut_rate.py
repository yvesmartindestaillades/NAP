from email.mime import base
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from os.path import exists, dirname
import os, sys
import numpy as np
import seaborn as sns
import json
from dreem_nap import data_wrangler, data_manip, database, plot, utils
from dreem_nap.study import Study


# Set your root folder for the database (at the moment, keep Yves)
folder = 'Yves'

path_to_data = 'data'

# Firebase credentials file
firebase_credentials_file = f"{path_to_data}/credentials_firebase.json"
with open(firebase_credentials_file) as file:
    firebase_credentials = json.load(file)

# Select your study
study_name = 'all samples' 

## Set your base coverage high-pass filter value
min_bases_cov = 1000 

# Set the resolution for the plots
mpl.rcParams['figure.dpi'] = 600 # the highest the resolution, the slowest the plotting
mpl.rcParams["figure.figsize"] = [25,7]
#plt.rcParams["figure.autolayout"] = True

# Depending on the study you select, you'll get a series of samples. You can also create new studies using this dictionary.
# Here's an example.
studies = data_wrangler.load_studies( f"{path_to_data}/samples.csv")

study = studies.loc[study_name]


## Pickle files to process and to push to Firebase
pickles_list = []# study.samples # Can be samples if you want to process the samples from your study

pickles = {pickle:  f"{path_to_data}/DREEM/{pickle}/mutation_histos.p" for pickle in pickles_list}

# Indicate the location of your RNA structure file
RNAstructureFile =  f"{path_to_data}/RNAstructureFile.csv"

# Default location for your local database (JSON file)
json_file =  f"{path_to_data}/db.json"

# If the user gives some new pickles files, push them to the firebase, then pull the entire firebase
if len(pickles):
    data_wrangler.push_samples_to_firebase(pickles = pickles,
                        RNAstructureFile = RNAstructureFile,
                        firebase_credentials = firebase_credentials,
                        min_bases_cov = min_bases_cov, 
                        folder=folder)

# Pull the firebase
#df_database = database.load(study=study, folder=folder)

#data_wrangler.dump_dict_json(JSONFileDict=json_file, df=df_database)
df_database = data_wrangler.json_load(json_file)

# Clean and reformat the dataset
df, df_full = data_wrangler.clean_dataset(df_database=df_database,
                                             study=study)
print(f"{df.groupby('construct')['samples_covered'].count().count()} constructs have more than {min_bases_cov} reads for each base of their ROI on each sample")
        


def study_base_wise_mut_rate(df:pd.DataFrame, study:Study, construct:int, bases = ['A','C'], scale_x = 'lin', figsize=(24,10))->None:
    """Generate line-plots of each base's mutation rate w.r.t a study's conditions, for a specific construct.

    Args:
        df (pd.DataFrame): dataframe of interest.
        study (Study): class containing relevant information about the series of sample that you want to use.
        construct (int): construct of interest.
        bases (list[str]): bases to display, sublist of ['A','C','G','T']
        scale_x (str): linear 'lin' or log 'log'
        figsize (Tuple(int,int)): size of the plotted figure.
    """

    df_paired, df_not_paired = pd.DataFrame(), pd.DataFrame()
    for samp in study.samples:
        df_roi = data_manip.get_roi_info(df, samp, construct, bases)
        df_paired = pd.concat((df_paired, 
                              df_roi['mut_rate'].xs(True, level='paired').reset_index().set_index('index')
                              .drop(columns=['base','roi_structure_comparison']).transpose()))
        df_not_paired = pd.concat((df_not_paired, 
                                   df_roi['mut_rate'].xs(False, level='paired').reset_index().set_index('index')
                                   .drop(columns=['base','roi_structure_comparison']).transpose()))

    df_paired, df_not_paired = df_paired.set_index(pd.Series(study.conditions)), df_not_paired.set_index(pd.Series(study.conditions))

    # Plot it
    fig = plt.figure()
    fig.suptitle(f"Construct {construct}, {study.name}", fontsize=16)
    ax1, ax2 = plt.subplot(121), plt.subplot(122)

    df_paired.plot(figsize=figsize,
                    logx={'lin':False, 'log':True}[scale_x],
                    ax=ax1, 
                    sharey=True, 
                    title='Paired bases',  
                    xlabel=f"{study.title}",
                    ylabel="Mutation rate")

    df_not_paired.plot(figsize=figsize,
                    logx={'lin':False, 'log':True}[scale_x],
                    ax=ax2, 
                    sharey=True, 
                    title='Unpaired bases', 
                    xlabel=f"{study.title}",
                    ylabel="Mutation rate")

  #  plt.tight_layout()
for stu in studies.iterrows():
    for construct in df.construct.unique():
        study=Study().from_dict(stu[1].to_dict())
        if study.name == 'all samples':
            continue
        study_base_wise_mut_rate(df=df,
                                study=study,
                                construct=construct)
        utils.save_fig(path= f"{path_to_data}/figs/date/Base-wise mutation along a study/{study.name}", 
                       title=f"{construct}_{study.name}.png")
        plt.close()