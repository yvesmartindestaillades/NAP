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
    sys.path.append(dirname('../libs/dreem/dreem')) 
except:
    "If dreem isn't installed on your computer, the code won't run"

from libs import dreem
from scipy.stats import linregress
from matplotlib.offsetbox import AnchoredText
from nap import utils, firebase

CONST_R = 1.98720425864083E-3 #Kcal.K^-1.mol^-1
CONST_T = 310.15 #KELVINS



def clean_dataset(df_rough, tubes, min_bases_cov):
    # Only keep desired pickle files
    df_full = df_rough[df_rough['tube'].isin(tubes)]

    # Check how many tubes reach 1000 reads on each base for a given construct
    df_full['tubes_covered'] = pd.Series(dtype=int)
    for construct in df_full.groupby('construct'):
        df_full['tubes_covered'].loc[construct[1].index] = construct[1]['full_sequence'].count()

    # Only keep constructs that reach 1000 reads in every tube    
    df = df_full[df_full['tubes_covered'] == len(tubes)].reset_index().drop(columns='index')

    number_of_void_dropped = (utils.added_content_per_construct(df, 'roi_deltaG' )=='void').apply(int).sum()
    print(f"{number_of_void_dropped} constructs were dropped because deltaG was 'void'")
    df = df[df['roi_deltaG'] != 'void']

    df = df.astype(dtype={'tube': str, 'construct':int, 'roi_sequence':str, 'full_sequence':str, 'roi_start_index':int,
    'roi_end_index':int, 'roi_deltaG':float, 'full_deltaG':float,
    'roi_structure_comparison':str, 'full_structure':str, 'data_type':str,
    'num_reads':int, 'num_aligned':int, 'num_of_mutations':object, 'mut_bases':object,
    'info_bases':object, 'del_bases':object, 'ins_bases':object, 'cov_bases':object, 'start':int, 'end':int,
    'mod_bases_A':object, 'mod_bases_C':object, 'mod_bases_G':object, 'mod_bases_T':object,
    'skips_low_mapq':int, 'skips_short_read':int, 'skips_too_many_muts':int,
    'cov_bases_roi':int, 'cov_bases_sec_half':int, 'tubes_covered':int,
    'sub-library':str, 'flank':str})

    print(f"{df.groupby('construct')['tubes_covered'].count().count()} constructs have more than {min_bases_cov} reads for each base of their ROI on each tube")

    return df, df_full


def pickle2dict(mhs, dropAttribute):
    localDict = {}
    for construct in mhs:
        localDict[construct] = mhs[construct].__dict__
        for attribute in dropAttribute:
            del localDict[construct][attribute]

        np_arrays = ['mut_bases', 'info_bases', 'del_bases', 'ins_bases',
                    'cov_bases']
        for array in np_arrays:
            localDict[construct][array] = tuple(localDict[construct][array])
        
        np_bases_arrays = ['A', 'C', 'G', 'T']
        for array in np_bases_arrays:
            localDict[construct]['mod_bases_'+array] = tuple(localDict[construct]['mod_bases'][array])
        del localDict[construct]['mod_bases']

        skips = ['low_mapq', 'short_read', 'too_many_muts']
        for sk in skips:
            localDict[construct]['skips_'+sk] = localDict[construct]['skips'][sk]
        del localDict[construct]['skips']
    return localDict


def generate_pickles(path_to_data,  pickles_list= None, letters_boundaries=['B','A'], number_boundaries=[1,0], remove_pickles=[]):
    list_of_pickles, pickles = pickles_list, {}
    alphabet = list(string.ascii_uppercase)
    for letter in alphabet[alphabet.index(letters_boundaries[0]):alphabet.index(letters_boundaries[1])+1]:
        for number in range(number_boundaries[0],number_boundaries[1]+1):
            list_of_pickles.append(letter+str(number))
    
    for items in remove_pickles:  
        try:
            list_of_pickles.remove(items)
        except:
            continue
            
    for pickle in list_of_pickles:
        pickles[pickle] = f"{path_to_data}/{pickle}/mutation_histos.p"

    return pickles

    
def dump_string_json(JSONFileString, df):
    print(f"Dumping df as a string to a JSON file {JSONFileString}")
    with open(JSONFileString, 'w') as outfile:
        json.dump(df.to_json(orient='index'), outfile) 
    print("Done!")

def load_string_json(JSONFileString):
    print("Load from JSON file")
    with open(JSONFileString) as json_file:
        my_json = json.load(json_file)
        df = pd.read_json(my_json, orient='index')
    print("Done!")
    return df

def dump_dict_json(JSONFileDict, df):
    print(f"Dumping df as a dict to a JSON file {JSONFileDict}")
    with open(JSONFileDict, 'w') as outfile:
        this_dict = df.set_index(['tube', 'construct']).groupby(level=0)\
            .apply(lambda d: d.reset_index().set_index('construct').to_dict(orient='index')).to_dict()
        json.dump(this_dict , outfile)
    print("Done!")
    return this_dict

def load_dict_json(JSONFileDict):
    print("Load from dict-type JSON file")
    with open(JSONFileDict) as json_file:
        dictionary = json.load(json_file)
    # dictionary = pd.DataFrame.from_dict(my_json, orient='columns')
        df = pd.DataFrame.from_dict({(i,j): dictionary[i][j] 
                                for i in dictionary.keys() 
                                for j in dictionary[i].keys()},
                            orient='columns').transpose()\
                            .drop(columns='tube')\
                            .dropna()\
                            .reset_index()\
                            .rename({'level_0': 'tube', 'level_1': 'construct'}, axis=1)
    print("Done!")
    return df


def push_pickles_to_firebase(pickles, RNAstructureFile, min_bases_cov, username, print_end=' '):
    # Load additional content
    df_additional_content = pd.read_csv(RNAstructureFile)
    df_additional_content.construct = df_additional_content.construct.astype(int).astype(str)
    df_additional_content.full_sequence = df_additional_content.full_sequence.apply(lambda seq: seq.replace('U','T'))

    print('Push pickles to firebase!')
    for count, tube in enumerate(pickles):
        # Load a tube from a pickle file
        mhs = pickle.load(open(pickles[tube], "rb"))

        df_tube = pd.DataFrame.from_dict(pickle2dict(mhs, dropAttribute = ['structure','_MutationHistogram__bases','sequence']),
                orient='index').rename(columns={'name':'construct'})

        # Merge with additional content (excel sheet content) and check the data sanity by sequences comparison
        df_temp = pd.merge(df_additional_content, df_tube, how='inner', on='construct').reset_index()
        #  assert not df_temp.apply(lambda row: not (str(row['full_sequence']).replace('U','T') in str(row['sequence'])) ,axis=1).sum(), "A sequence didn't match in the fusion"          
        df_temp = df_temp.drop(columns='index')

        # Count base coverage in the ROI and in the second half            
        df_temp['cov_bases_roi'] = df_temp.apply(lambda row: np.array(np.array(row['cov_bases'])[int(row['roi_start_index']):int(row['roi_end_index'])]).min(), axis=1)
        df_temp['cov_bases_sec_half'] = df_temp.apply(lambda row: np.array(np.array(row['cov_bases'])[int(len(row['cov_bases'])/2):]).min(), axis=1)

        # Filter out the constructs that don't reach 1000 reads for each base of the ROI 
        df_temp = df_temp[df_temp['cov_bases_roi'] >= min_bases_cov]
        
        df_temp = df_temp.astype(dtype={'construct':int, 'roi_sequence':str, 'full_sequence':str, 'roi_start_index':int,
        'roi_end_index':int, 'roi_structure_comparison':str, 'full_structure':str, 'data_type':str,
        'num_reads':int, 'num_aligned':int, 'num_of_mutations':object, 'mut_bases':object,
        'info_bases':object, 'del_bases':object, 'ins_bases':object, 'cov_bases':object, 'start':int, 'end':int,
        'mod_bases_A':object, 'mod_bases_C':object, 'mod_bases_G':object, 'mod_bases_T':object,
        'skips_low_mapq':int, 'skips_short_read':int, 'skips_too_many_muts':int,
        'cov_bases_roi':int, 'cov_bases_sec_half':int, 'sub-library':str, 'flank':str})

        df_temp = df_temp.set_index('construct')

        # Push this tube to firebase
        firebase.push(df_temp.to_dict(orient='index'), ref=tube, username=username, verbose= not bool(count))

        # Give yourself hope to wait by showing the progress
        print(tube, end=print_end)
    print('Done!')