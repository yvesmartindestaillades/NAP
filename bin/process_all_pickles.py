import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from os.path import exists, dirname
import os, sys

script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..' )
sys.path.append( mymodule_dir )

from NAP.nap import *
try:
    sys.path.append(dirname('libs/dreem/dreem')) 
except:
    "If dreem isn't installed on your computer, the code won't run"
    
# Set your username for the database (at the moment, keep Yves)
username = 'Gabe'

## Set your base coverage high-pass filter value
min_bases_cov = 1000 

## Pickle files to process and to push to Firebase
pickles_list =  [ele for ele in [f"{a}{b}" for a in string.ascii_uppercase[0:8] for b in range(1,11)] \
                 if ele not in ['C3','C10','D10','E10','F10','G10','H10', 'E4']]\
                 + ['C5_realignment_v3']

pickles = data_wrangler.generate_pickles(path_to_data='data/FULLSET',
                                         pickles_list=pickles_list)

# Indicate the location of your RNA structure file
RNAstructureFile = 'data/RNAstructureFile.csv'

# Default location for your local database (JSON file)
json_file = 'data/db.json'

# If the user gives some new pickles files, push them to the firebase, then pull the entire firebase
if len(pickles):
    firebase.push_pickles(pickles = pickles,
                        RNAstructureFile = RNAstructureFile,
                        min_bases_cov = min_bases_cov, 
                        username=username, 
                        print_end='\n')