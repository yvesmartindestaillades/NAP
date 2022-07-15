import pandas as pd
import pickle
import json
import numpy as np
from os.path import dirname
import sys, os
from typing import Tuple, List

from dreem_nap import database, data_manip
from dreem_nap.study import Study
from dreem.bit_vector import MutationHistogram


def json_dump(df:pd.DataFrame, json_file:str, verbose:bool=True)->None:
    """A simple function to dump a Pandas dataframe into a (sample, construct)-wise indexed json file.

    Args:
        df (pd.DataFrame): the Pandas dataframe that you want to dump into a json.
        json_file (str): a string containing the relative path + name of the json to dump.
        verbose (bool, optional):  print relevant information. Defaults to True.
    Returns:
        None
    """

    if verbose: print(f"Dumping df as a dict to a JSON file {json_file}")
    with open(json_file, 'w') as outfile:
        this_dict = df.set_index(['samp', 'construct']).groupby(level=0)\
            .apply(lambda d: d.reset_index().set_index('construct').to_dict(orient='index')).to_dict()
        json.dump(this_dict , outfile)
    if verbose: print("Done!")




def json_load(json_file:str, verbose:bool = True)->pd.DataFrame:
    """A simple function to load data from a (sample, construct)-wise indexed json file.
    
    Args:
        json_file: a string containing the relative path + name of the json to load.
        verbose: print relevant information

    Returns:
        Loaded Pandas dataframe.
    """
    if verbose: print("Load from dict-type JSON file")
    with open(json_file) as file:
        dictionary = json.load(file)
        df = pd.DataFrame.from_dict({(i,j): dictionary[i][j] 
                                for i in dictionary.keys() 
                                for j in dictionary[i].keys()},
                            orient='columns').transpose()\
                            .drop(columns='samp')\
                            .dropna()\
                            .reset_index()\
                            .rename({'level_0': 'samp', 'level_1': 'construct'}, axis=1)
    if verbose: print("Done!")
    return df


def mhs2dict(mhs:MutationHistogram, drop_attribute:List[str])->dict:
    """Turns the output of DREEM into a construct-wise index dictionary.

    Args:
        mhs (_MutationHistogram): one sample's content under DREEM's MutationHistogram class format. 
        drop_attribute (List[str]): a list of attributes from MutationHistogram class that you don't want into your dictionary
    
    Returns:
        A dictionary form of the MutationHistogram class.
    """
    sample_dict = {}
    for construct in mhs:
        sample_dict[construct] = mhs[construct].__dict__
        for attribute in drop_attribute:
            del sample_dict[construct][attribute]

        np_arrays = ['mut_bases', 'info_bases', 'del_bases', 'ins_bases',
                    'cov_bases']
        for array in np_arrays:
            sample_dict[construct][array] = tuple(sample_dict[construct][array])
        
        np_bases_arrays = ['A', 'C', 'G', 'T']
        for array in np_bases_arrays:
            sample_dict[construct]['mod_bases_'+array] = tuple(sample_dict[construct]['mod_bases'][array])
        del sample_dict[construct]['mod_bases']

        skips = ['low_mapq', 'short_read', 'too_many_muts']
        for sk in skips:
            sample_dict[construct]['skips_'+sk] = sample_dict[construct]['skips'][sk]
        del sample_dict[construct]['skips']
    return sample_dict


def push_samples_to_firebase(pickles:dict, RNAstructureFile:str, min_bases_cov:int, firebase_credentials:dict, folder:str, verbose:bool=True, print_end:str=' ')->None:
    """Pushes new samples to Firebase.

    DREEM module outputs MutationHistogram objects, compressed under the pickle format. 
    For each pickle in pickles, this function turns it into dictionaries, filters out unvalid constructs and pushes the resulting (sample, construct)-wise indexed dictionary to the Firebase.
    The construct high-pass filter filters out a construct if which at least one base in its Region of Interest (ROI) doesn't reach `min_bases_cov` of base coverage.
    
    Args:
        pickles (dict): a dictionary with names of the samples and location of their pickle file: {'name of the sample': 'path to the file'}
        RNAstructureFile (str): string containing the name of a csv file with additional content. Its columns are: ['construct', 'roi_sequence', 'full_sequence', 'roi_start_index', 'roi_end_index', 'roi_deltaG', 'full_deltaG', 'roi_structure_comparison', 'full_structure', 'flank', 'sub-library']
        firebase_credentials (dict): Firebase credentials to access the database.
        min_bases_cov (int): int type. Mutation rates of bases below this threshold will be considered irrelevant.
        folder (str): where to push the data in the firebase.
        verbose (bool): print relevant information

    Returns:
        None
    """
    
    # Load additional content
    df_additional_content = pd.read_csv(RNAstructureFile)
    df_additional_content.construct = df_additional_content.construct.astype(int).astype(str)
    df_additional_content.full_sequence = df_additional_content.full_sequence.apply(lambda seq: seq.replace('U','T'))

    if verbose: print('Push pickles to firebase!')
    for count, samp in enumerate(pickles):
        # Load a sample from a pickle file
        mhs = pickle.load(open(pickles[samp], "rb"))

        df_sample = pd.DataFrame.from_dict(mhs2dict(mhs, drop_attribute = ['structure','_MutationHistogram__bases','sequence']),
                orient='index').rename(columns={'name':'construct'})

        # Merge with additional content (excel sheet content) and check the data sanity by sequences comparison
        df_temp = pd.merge(df_additional_content, df_sample, how='inner', on='construct').reset_index()
        #  assert not df_temp.apply(lambda row: not (str(row['full_sequence']).replace('U','T') in str(row['sequence'])) ,axis=1).sum(), "A sequence didn't match in the fusion"          
        df_temp = df_temp.drop(columns='index')

        # Count base coverage in the ROI and in the second half            
        df_temp['cov_bases_roi'] = df_temp.apply(lambda row: np.array(np.array(row['cov_bases'])[int(row['roi_start_index']):int(row['roi_end_index'])]).min(), axis=1)
        df_temp['cov_bases_sec_half'] = df_temp.apply(lambda row: np.array(np.array(row['cov_bases'])[int(len(row['cov_bases'])/2):]).min(), axis=1)

        # Filter out the constructs that don't reach 1000 reads for each base of the ROI 
        df_temp = df_temp[df_temp['cov_bases_roi'] >= min_bases_cov]
        df_temp['min_bases_cov'] = min_bases_cov

        df_temp['base_pairing_prob'] = df_temp['base_pairing_prob'].apply(lambda constr: [float(p) for p in constr.split(" ")].insert(0,0))

        df_temp = df_temp.astype(dtype={'construct':int, 'roi_sequence':str, 'full_sequence':str, 'roi_start_index':int,
        'roi_end_index':int, 'roi_structure_comparison':str, 'full_structure':str, 'data_type':str,
        'num_reads':int, 'num_aligned':int, 'num_of_mutations':object, 'mut_bases':object,
        'info_bases':object, 'del_bases':object, 'ins_bases':object, 'cov_bases':object, 'start':int, 'end':int,
        'mod_bases_A':object, 'mod_bases_C':object, 'mod_bases_G':object, 'mod_bases_T':object,
        'skips_low_mapq':int, 'skips_short_read':int, 'skips_too_many_muts':int,
        'cov_bases_roi':int, 'cov_bases_sec_half':int, 'sub-library':str, 'flank':str, 'min_bases_cov':int})

        df_temp = df_temp.set_index('construct')

        # Push this sample to firebase
        database.connect(firebase_credentials=firebase_credentials, verbose=False)
        database.push(df_temp.to_dict(orient='index'), ref=samp, folder=folder)

        # Give yourself hope to wait by showing the progress
        print(samp, end=print_end)
    print('Done!')



def filter_constructs_study_wise(df_database:pd.DataFrame, studies, verbose:bool = True):
    """ Removes constructs that aren't present for each sample of a study. 

    Args:
        df_database (pd.Dataframe): the dataframe downloaded from the Firebase.
        studies (dict[Study]): the studies you want to filter your constructs with. Format is {name:Study} 
        verbose: print relevant information

    Returns:
        df (pd.DataFrame): a dataframe with deltaG != void
        df_void (pd.DataFrame): a dataframe with deltaG == void
        df_full (pd.DataFrame): the merge of df and df_void

    """

    # Only keep desired pickle files
    samples = []
    for s in pd.DataFrame.from_dict({study: studies[study].to_dict()  for study in studies}, orient='index')['samples']:
        samples += s
    df_full = df_database[df_database['samp'].isin(samples)]

    # Only keep desired pickle files
    df_full = df_database[df_database['samp'].isin(samples)]

    # Check how many samples reach 1000 reads on each base for a given construct
    df_full['samples_covered'] = pd.Series(dtype=int)
    for construct in df_full.groupby('construct'):
        df_full['samples_covered'].loc[construct[1].index] = construct[1]['full_sequence'].count()

    # Only keep constructs that reach 1000 reads in every sample    
    # Check how many samples reach 1000 reads on each base for a given construct
    
    df = pd.DataFrame()

    for study in studies:
        cont = 0
        for doubles in ['replicate','all samples', 'TO DO']:
            if studies[study].name in doubles:
                cont = 1
        if cont:
            continue
        df_loc = df_full[df_full['samp'].isin(studies[study].samples)].set_index('construct')
        # print(df_loc.groupby('construct').count()['samp'])
        df_loc = df_loc[df_loc.groupby('construct').count()['samp'] == len(studies[study].samples)]
        df = pd.concat((df, df_loc))
    df = df.reset_index() 
    df = df.loc[df.astype(str).drop_duplicates().index]


    number_of_void_dropped = (data_manip.get_construct_attribute(df, 'roi_deltaG' )=='void').apply(int).sum()
    if verbose: print(f"{number_of_void_dropped} constructs were dropped because deltaG was 'void'")
    df_full = df.copy()
    
    df = df[df['roi_deltaG'] != 'void']
    df_void = df[df['roi_deltaG'] == 'void'].copy()



    df = df.astype(dtype={'samp': str, 'construct':int, 'roi_sequence':str, 'full_sequence':str, 'roi_start_index':int,
    'roi_end_index':int, 'roi_deltaG':float, 'full_deltaG':float,
    'roi_structure_comparison':str, 'full_structure':str, 'data_type':str,
    'num_reads':int, 'num_aligned':int, 'num_of_mutations':object, 'mut_bases':object,
    'info_bases':object, 'del_bases':object, 'ins_bases':object, 'cov_bases':object, 'start':int, 'end':int,
    'mod_bases_A':object, 'mod_bases_C':object, 'mod_bases_G':object, 'mod_bases_T':object,
    'skips_low_mapq':int, 'skips_short_read':int, 'skips_too_many_muts':int,
    'cov_bases_roi':int, 'cov_bases_sec_half':int, 'samples_covered':int,
    'sub-library':str, 'flank':str})

    return df, df_void, df_full





def load_studies(studies_file_path:str)->pd.DataFrame:
    """Read formatted file with samples, and turn it into a dataframe containing studies.

    Args:
        studies_file_path (str): path+title of the csv file containing the samples.

    Returns:
        (pd.DataFrame): studies of the csv file, indexed by study.

    Example:
        >>> from dreem_nap import data_wrangler     
        >>> study_file = 'samples.csv'
        >>> samples = ['"name","description","samples","conditions","title"',
        ... ',,,,',
        ... '"salt","Change the Na concentration","A6",0.15,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","B6",0.3,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","C6",0.6,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","D6",1,"Na quantity [M]"',
        ... '"salt","Change the Na concentration","E6",1.2,"Na quantity [M]"',
        ... ',,,,',
        ... '"spermidine","Change the Spermidine concentration","B3",0.01,"Spermidine quantity [mM]"',
        ... '"spermidine","Change the Spermidine concentration","D3",1,"Spermidine quantity [mM]"',
        ... '"spermidine","Change the Spermidine concentration","E3",10,"Spermidine quantity [mM]"',
        ... '"spermidine","Change the Spermidine concentration","F3",100,"Spermidine quantity [mM]"']
        >>> with open(study_file, 'w') as f:
        ...     f.writelines(f"{l}\\n" for l in samples)
        >>> # Your studies under the shape of a dataframe
        >>> df_studies = data_wrangler.load_studies(study_file) 
        >>> df_studies
                                name                           description                     samples                          title                        conditions
        salt                    salt           Change the Na concentration        [A6, B6, C6, D6, E6]                Na quantity [M]        [0.15, 0.3, 0.6, 1.0, 1.2]
        spermidine        spermidine   Change the Spermidine concentration            [B3, D3, E3, F3]       Spermidine quantity [mM]          [0.01, 1.0, 10.0, 100.0]
        >>> temp = df_studies.to_dict(orient='index')
        >>> # Your studies under the shape of a dictionary of Study
        >>> studies = {study: Study().from_dict(temp[study])  for study in temp}
        >>> print(f"Here are the studies: {studies.keys()}")
        dict_keys(['salt', 'spermidine'])
        >>> study_name = 'salt' 
        >>> study = studies[study_name] 
        >>> print(f"Here is your study {study.to_dict()}" )
        Here is your study {'name': 'salt', 'description': 'Change the Na concentration', 'samples': ['A6', 'B6', 'C6', 'D6', 'E6'], 'title': 'Na quantity [M]', 'conditions': [0.15, 0.3, 0.6, 1.0, 1.2]}
    """

    studies_dict, studies_data = {}, pd.read_csv(studies_file_path)

    for col in studies_data.groupby('name')[Study().attr_list]:
        solo_item = lambda x: x[0] if len(set(x)) == 1 else x  
        studies_dict[col[0]] = {attr: solo_item(list(col[1][attr])) for attr in (Study().attr_list)} 
    
    return pd.DataFrame.from_dict(studies_dict, orient='index')
