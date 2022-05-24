# %%
import pickle
import json
import pandas as pd
import plotly.express as px
import numpy as np
import plotly.graph_objects as go

def progBarIter(progBar, counter, character):
    new = list(progBar)
    new[counter] = character
    return ''.join(new)

class mutationHistogramManipulation():

    def pickle2dict(mhs, dropAttribute):
        localDict = {}
        for construct in mhs:
            localDict[construct] = mhs[construct].__dict__
            for attribute in dropAttribute:
                del localDict[construct][attribute]

            # Turn mutationHistogram numpy arrays into lists
            # If the bit_vector class is updated in the dreem pipeline, 
            # this code could be replaced by:
            # localDict[construct] = mhs[construct].np_arrays_to_lists.__dict__
            # See bit_vector.py

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

    def fusionPickles(pickles, dropAttribute):
        out_pickles = pickles.copy()
        globalDict, progBar = {}, '_'*len(pickles)
        for counter, pickleFile in enumerate(pickles.keys()):
            try:
                mhs = pickle.load(open(pickles[pickleFile], "rb"))
            except:
                out_pickles.pop(pickleFile)
                progBar = progBarIter(progBar, counter, pickleFile )
                print(progBar)
                continue

            localDict = mutationHistogramManipulation.pickle2dict(mhs, dropAttribute)

            progBar = progBarIter(progBar, counter, '#')
            print(progBar)
            globalDict[pickleFile] = localDict
        return globalDict, out_pickles

    def dictToJSON(dict, fileName):
        print("Dumping the dict to a JSON file...")
        with open(fileName, 'w') as outfile:
            json.dump(dict, outfile) 
        print("Done!") 

    def JSONtoDict(JSONFile):
        print("Load from JSON file")
        with open(JSONFile) as json_file:
            data = json.load(json_file)
        print("Done!")
        return data

    def aggregateTwoDict(main_dict, aggregate_dict, aggregate_columns):
        for tube in main_dict.keys():
            for construct in main_dict[tube].keys(): 
                # Check that the sequence is good
                assert aggregate_dict['full_sequence'][construct].replace('U','T') in main_dict[tube][construct]['sequence'],\
                        print(f"Failing index is {construct}")
                # Give new attributes
                for attribute in aggregate_columns:
                    main_dict[tube][construct][attribute] = aggregate_dict[attribute][construct]
        return main_dict
# %%
