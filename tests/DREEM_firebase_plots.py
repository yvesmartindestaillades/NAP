import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from os.path import exists, dirname
import os, sys
from numpy import random
import seaborn as sns 
import shutil

#sys.path.append(os.path.abspath(".."))

sys.path.append(os.path.abspath(""))
from NAP import data_wrangler, utils, firebase, plot

try:
    shutil.rmtree('data/test_output')
except:
    'no previous test_output dir'
os.mkdir('data/test_output')

# Set your username for the database (at the moment, keep Yves)
username = 'test'

# Select your study
study = 'test' 

## Set your base coverage high-pass filter value
min_bases_cov = 1000 

# Set the resolution for the plots
mpl.rcParams['figure.dpi'] = 50 # the highest the resolution, the slowest the plotting
plt.rcParams["figure.figsize"] = [25, 7] # you still want regular-sized plots

# Depending on the study you select, you'll get a series of tubes. You can also create new studies using this dictionary.
# Here's an example.
tubes_per_study = {   
    'test':             ['A6', 'D6','F3','G1'],
    }

tubes = tubes_per_study[study]



## Pickle files to process and to push to Firebase
# Can be tubes if you want to process the tubes from your study, or [] if they are already on the database 
pickles_list = tubes

pickles = data_wrangler.generate_pickles(path_to_data='data/FULLSET',
                                         pickles_list=pickles_list)

# Indicate the location of your RNA structure file
RNAstructureFile = 'data/RNAstructureFile.csv'

# If the user gives some new pickles files, push them to the firebase, then pull the entire firebase
if len(pickles): 
    data_wrangler.push_pickles_to_firebase(pickles = pickles,
                                            RNAstructureFile = RNAstructureFile,
                                            min_bases_cov = min_bases_cov, 
                                            username=username)


# Pull the firebase
df_rough = firebase.load(tubes=tubes, username=username)

assert len(df_rough), "Couldn't pull from Firebase"

# Clean and reformat the dataset
df, df_full = data_wrangler.clean_dataset(df_rough=df_rough,
                                             tubes=tubes, 
                                             min_bases_cov=min_bases_cov)

assert len(df.tube.unique()) == len(tubes), "Not all tubes were loaded"

assert len(df.construct.unique()) > 10, "Less than 10 valid constructs in the dataset"

# Test data quality plots
plot.base_coverage_for_all_constructs(df=df_full, 
                                      min_bases_cov=min_bases_cov)
plot.save_fig(path='data/test_output',
                title='base_coverage_for_all_constructs')
plt.close()

plot.random_9_base_coverage(df=df, 
                            min_bases_cov=min_bases_cov)
plot.save_fig(path='data/test_output',
                title='random_9_base_coverage')
plt.close()

tube, construct = utils.rand_tube_construct(df)
plot.base_coverage(df, tube, construct , min_bases_cov=min_bases_cov)
plot.save_fig(path='data/test_output',
                title='base_coverage')
plt.close()

plot.heatmap(df = df, 
             column="cov_bases_roi")
plot.save_fig(path='data/test_output',
                title='heatmap')
plt.close()

# Test data analysis plots
tube, construct = utils.rand_tube_construct(df)
plot.mutation_rate(df=df,
                    tube=tube,
                    construct=construct,
                    plot_type='sequence',
                    index='index')
plot.save_fig(path=f"data/test_output", 
            title=f"mutation_rate_{tube}_{construct}")
plt.close()

tube, _ = utils.rand_tube_construct(df)
plot.deltaG(df=df, tube=tube)
plot.save_fig(path=f"data/test_output", 
            title=f"deltaG_{tube}")
plt.close()

_ , construct = utils.rand_tube_construct(df)
df_global_corr = plot.correlation_n_tubes(df, tubes, construct)
plot.save_fig(path=f"data/test_output", 
                title=f"correlation_{study}_{construct}")
plt.title(f"Correlation_for_{study}")
plt.close()

_ , constructs = utils.rand_tube_construct(df, n_constructs=3)
for plt_type in ['r_value', 'slope']:
    pivot = df_global_corr.pivot("tube_0","tube_1", plt_type).astype(float)
    f, ax = plt.subplots(figsize=(28, 10))
    sns.heatmap(pivot, annot=False, linewidths=0, ax=ax)#, norm=LinNorm())
    plt.title(f"{plt_type} of the correlation between tubes for study: {study}, constructs: {constructs}")
    plot.save_fig(path=f"data/test_output", 
                    title=f"correlation_{plt_type}_{study}_{constructs}")
    plt.close()

utils.columns_to_csv(df=df,
                   tubes=tubes,
                   columns=['tube', 'construct','full_sequence','roi_sequence','mut_bases','info_bases'],
                   title=f"seq_and_reactivity_{study}",
                   path='data/test_output')

utils.deltaG_vs_construct_to_csv(df=df,    
                                 title=f"deltaG_vs_construct.csv", 
                                 path = f"data/test_output", 
                                 tubes=tubes)

