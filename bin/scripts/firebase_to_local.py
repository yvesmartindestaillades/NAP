import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from os.path import exists, dirname
import os, sys
sys.path.append(os.path.abspath(""))

try:
    sys.path.append(dirname('libs/dreem/dreem')) 
except:
    "If dreem isn't installed on your computer, the code won't run"
    
# Set your username for the database (at the moment, keep Yves)
username = 'Yves'

## Set your base coverage high-pass filter value
min_bases_cov = 1000 

## Pickle files to process and to pull from Firebase
tubes =  [ele for ele in [f"{a}{b}" for a in string.ascii_uppercase[0:8] for b in range(1,11)] \
                 if ele not in ['C3','C10','D10','E10','F10','G10','H10', 'E4']]\
                 + ['C5_realignment_v3']

# Default location for your local database (JSON file)
json_file = 'data/db.json'

# Pull the firebase
df_rough = firebase.load(tubes=tubes, username=username)

# Clean and reformat the dataset
df, df_full = data_wrangler.clean_dataset(df_rough=df_rough,
                                             tubes=tubes, 
                                             min_bases_cov=min_bases_cov)

# Dump to JSON 
data_wrangler.dump_dict_json(json_file, df_full)