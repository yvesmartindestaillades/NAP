import pandas as pd
import numpy as np

def __find_base_in_sequence(sequence, base_type):
    return [i for i, base in enumerate(sequence) if base in base_type]

def __filtered_index(row, base_index, base_type, base_pairing, RNAstructure_use_DMS, RNAstructure_use_temp):
    index = list(range(len(row['sequence'])))
    if base_index is not None:
        if isinstance(base_index, int):
            index = [base_index]
        if isinstance(base_index, list):
            index = base_index
        if isinstance(base_index, str):
            assert (count:=row['sequence'].count(base_index)) == 1, f"{count} sequences {base_index} found in sequence of sample {row['sample']} construct {row['construct']} section {row['section']} cluster {row['cluster']} (sequence: {row['sequence']})"
            temp_idx = row['sequence'].find(base_index)
            index = list(range(temp_idx, temp_idx+len(base_index)))

    if base_type is not ['A','C','G','T']:
        index = list(set(index) & set(__find_base_in_sequence(row['sequence'], base_type)))
    
    if base_pairing is not None:
        base_pairing_switch = lambda idx, base, pairing: idx if (base=='.' and pairing) or (base in ['(',')'] and not pairing) else None 
        structure = [base_pairing_switch(i,base,base_pairing) for i, base in enumerate(row['structure'+('_DMS' if RNAstructure_use_DMS else '')+('_T' if RNAstructure_use_temp else '')])]
        index = list(set(index) & set([i for i in structure if i is not None]))
    return index

def get_df(df, sample=None, construct=None, section=None, cluster=None, min_cov_bases=0, base_index=None, base_type=['A','C','G','T'], base_pairing=None, RNAstructure_use_DMS=False, RNAstructure_use_temp=False, **kwargs)->pd.DataFrame:
    """Get a dataframe with filtered data

    Args:
        df (pd.Dataframe): Dataframe to filter
        sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
        construct (list, int, str, optional): Filter rows by construct (list of constructs or just a construct). Defaults to None.
        section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
        cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
        min_cov_bases (int, optional): Filter rows by a minimum threshold for base coverage. Defaults to 0.
        base_index (list, int, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
        base_type (list, str, optional): Filter per-base attributes (mut_rates, sequence, etc) by base type. Defaults to ['A','C','G','T'].
        base_pairing (bool, optional): Filter per-base attributes (mut_rates, sequence, etc) by predicted base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
        RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing. Defaults to False.
        RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing. Defaults to False.

    Returns:
        pd.Dataframe: a filtered dataframe according to the given parameters
    """


    df = df.copy()
    # filter mutation profiles
    df = df[df.worst_cov_bases >= min_cov_bases]
    for key, value in kwargs.items():
        locals()[key] = value
    mp_attr = ['sample', 'construct', 'section', 'cluster'] + list(kwargs.keys())

    for attr in mp_attr:
        if eval(attr) is not None:
            if (isinstance(eval(attr), list) or isinstance(eval(attr), tuple)):
                df = df[df[attr].isin(eval(attr))]
            else:
                df = df[df[attr] == eval(attr)]

    # filter base profiles
    df['filtered_index'] = df.apply(lambda row: __filtered_index(row, base_index, base_type, base_pairing, RNAstructure_use_DMS, RNAstructure_use_temp), axis=1)
    df = df[df.filtered_index.apply(lambda x: len(x) > 0)]

    bp_attr = ['sequence', 'mut_bases', 'info_bases','del_bases','ins_bases','cov_bases','mut_rates'] + \
        [c for c in df.columns.tolist() if (c.startswith('structure') or c.startswith('mod_bases') or c.startswith('poisson'))]
    for idx, row in df.iterrows():
        for attr in bp_attr:
            filtered_cell = [row[attr][i] for i in df.at[idx, 'filtered_index']]
            if type(row[attr]) == str:
                df.at[idx, attr] = ''.join(filtered_cell)
            else:
                df.at[idx, attr] = np.array(filtered_cell)
    return df

