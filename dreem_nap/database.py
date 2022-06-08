import pandas as pd
import firebase_admin
from firebase_admin import credentials
from firebase_admin import db
from os.path import exists, dirname
import os, sys
sys.path.append(os.path.abspath(""))
from dreem_nap import *
from typing import Tuple, List
from dreem_nap.study import Study

def connect(firebase_credentials:dict, verbose:bool = True)->None:
    """Initiate connection with the database.

    Args:
        credentials (dict): Firebase credentials to access the database.
        verbose (bool): print when connection was already set.
    """

    cred = credentials.Certificate(firebase_credentials)
    # Initialize the app with a service account, granting admin privileges
    # As an admin, the app has access to read and write all data, regradless of Security Rules

    try:
        default_app = firebase_admin.initialize_app(cred, {
            'databaseURL':'https://dreem-542b7-default-rtdb.firebaseio.com/'
            })
    except:
        if verbose: print('Re-used the previous Firebase connection')

def push(dict_df:pd.DataFrame, folder:str, ref:str)->None:
    """Push a dictionary to the database, in a given folder, for a given reference.
    
    Args:
        dict_df: the dictionary you want to push
        folder: string, root folder in the database. Corresponds to a user, a version, a project, etc.
        ref: string, path to w
    """

    ref_obj = db.reference(f"{folder}/{ref}")
    ref_obj.set(dict_df)


def load(folder:str, study:Study, verbose:bool = True)->pd.DataFrame:
    """Download a Pandas dataframe from the database.
    
    Args:
        folder: string, root folder in the database. Corresponds to a user, a version, a project, etc.
        study: class containing a list of samples that you want to use.
        verbose: print relevant information.
    
    Returns:
        Dataframe of the targeted samples.
    """
    samples = study.samples
    
    if verbose: print('Load data from database')
    connect()
    df = {}
    missed_samples = []
    for samp in samples:
        try:
            ref = db.reference(f"{folder}/{samp}")
            df[samp] = pd.DataFrame.from_dict(ref.get('/')[0], orient='index')
            if verbose: print(samp, end=' ')
        except:
            if verbose: print(f"\nsample {samp} not found on database")
            missed_samples.append(samp)

    if missed_samples != []:
        if verbose: print(f"samples {missed_samples} couldn't be loaded from database")

    df = pd.concat(df)
    df = df.reset_index().rename(columns={'level_0':'samp', 'level_1':'construct'})
    if verbose: print('Done!')
    return df