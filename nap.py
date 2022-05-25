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
from libs import dreem

CONST_R = 1.98720425864083E-3 #Kcal.K^-1.mol^-1
CONST_T = 310.15 #KELVINS


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

class turner_overthrow:
    def clean_dataset(df_rough, tubes, min_bases_cov):
        # Only keep desired pickle files
        df_full = df_rough[df_rough['tube'].isin(tubes)]

        # Check how many tubes reach 1000 reads on each base for a given construct
        df_full['tubes_covered'] = pd.Series(dtype=int)
        for construct in df_full.groupby('construct'):
            df_full['tubes_covered'].loc[construct[1].index] = construct[1]['sequence'].count()

        # Only keep constructs that reach 1000 reads in every tube    
        df = df_full[df_full['tubes_covered'] == len(tubes)].reset_index().drop(columns='index')

        number_of_void_dropped = (utils.added_content_per_construct(df, 'var_deltaG' )=='void').apply(int).sum()
        print(f"{number_of_void_dropped} constructs were dropped because deltaG was 'void'")
        df = df[df['var_deltaG'] != 'void']

        df = df.astype(dtype={'tube': str, 'construct':int, 'var_sequence':str, 'full_sequence':str, 'var_start_index':int,
        'var_end_index':int, 'var_deltaG':float, 'full_deltaG':float,
        'var_structure_comparison':str, 'full_structure':str, 'sequence':str, 'data_type':str,
        'num_reads':int, 'num_aligned':int, 'num_of_mutations':object, 'mut_bases':object,
        'info_bases':object, 'del_bases':object, 'ins_bases':object, 'cov_bases':object, 'start':int, 'end':int,
        'mod_bases_A':object, 'mod_bases_C':object, 'mod_bases_G':object, 'mod_bases_T':object,
        'skips_low_mapq':int, 'skips_short_read':int, 'skips_too_many_muts':int,
        'cov_bases_var':object, 'cov_bases_sec_half':object, 'tubes_covered':int,'var_structure_comparison':str})

        print(f"{df.groupby('construct')['tubes_covered'].count().count()} constructs have more than {min_bases_cov} reads for each base of their var part on each tube")

        return df, df_full


class data_wrangler:

    def generate_pickles(path_to_data,  pickles_list= None, letters_boundaries=[1,0], number_boundaries=['B','A'], remove_pickles=None):
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

    def push_pickles_to_firebase(pickles, RNAstructureFile, min_bases_cov, user):
        # Load additional content
        df_additional_content = pd.read_csv(RNAstructureFile)
        df_additional_content.construct = df_additional_content.construct.astype(str)
        df_additional_content = df_additional_content.set_index('construct')

        print('Push pickles to firebase!')
        for tube in pickles:
            # Load a tube from a pickle file
            mhs = pickle.load(open(pickles[tube], "rb"))

            df_tube = pd.DataFrame.from_dict(mutationHistogramManipulation.pickle2dict(mhs, dropAttribute = ['structure','_MutationHistogram__bases']),
                    orient='index').rename(columns={'name':'construct'}).set_index('construct')

            # Merge with additional content (excel sheet content) and check the data sanity by sequences comparison
            df_temp = pd.merge(df_additional_content, df_tube, left_index=True, right_index=True, how='inner')
            assert not df_temp.apply(lambda row: not (str(row['full_sequence']).replace('U','T') in str(row['sequence'])) ,axis=1).sum(), "A sequence didn't match in the fusion"

            # Count base coverage in the variable part and in the second half
            df_temp['cov_bases_var'] = df_temp.apply(lambda row: np.array(np.array(row['cov_bases'])[row['var_start_index']:row['var_end_index']]).min(),axis=1)
            df_temp['cov_bases_sec_half'] = df_temp.apply(lambda row: np.array(np.array(row['cov_bases'])[int(len(row['cov_bases'])/2):]).min(),axis=1)
            df_temp.head()

            # Filter out the constructs that don't reach 1000 reads for each base of the ROI 
            df_temp = df_temp[df_temp['cov_bases_var'] >= min_bases_cov]

            # Push this tube to firebase
            data_wrangler.push_data_to_firebase(df_temp.to_dict(orient='index'), ref=tube, user=user, verbose=False)

            # Give yourself hope to wait by showing the progress
            print(tube, end=' ')
        print('Done!')


    def push_data_to_firebase(dict_df, ref, user, verbose = True):
        if verbose: print("Push data to firebase")
        utils.connect_to_firebase()
        ref_obj = db.reference(f"{user}/{ref}")
        ref_obj.set(dict_df)
        if verbose: print("Done!")


    def load_data_from_firebase(tubes, user):
        print('No local file found, load data from Firebase')
        utils.connect_to_firebase()
        df = {}
        missed_tubes = []
        for tube in tubes:
            try:
                ref = db.reference(f"{user}/{tube}")
                df[tube] = pd.DataFrame.from_dict(ref.get('/')[0], orient='index')
                print(tube, end=' ')
            except:
                print(f"\nTube {tube} not found on Firebase")
                missed_tubes.append(tube)

        if missed_tubes != []:
            print(f"Tubes {missed_tubes} couldn't be loaded from Firebase")

        df = pd.concat(df)
        df = df.reset_index().rename(columns={'level_0':'tube', 'level_1':'construct'})
        print('Done!')
        return df
        
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

class plot:
    def tube_coverage_distribution(df):
        plt.figure(figsize=(25, 8))
        plt.plot(np.array(df['tubes_covered'].sort_values(ascending=False)))
        plt.xlabel('Constructs (sorted)')
        plt.ylabel('Amount of tubes covered')
        plt.grid()

    def valid_construct_per_tube(df, min_bases_cov, figsize=(25,8)):
        df.groupby('tube').count().reset_index().plot(kind='bar',x='tube', y='construct', figsize=figsize )
        plt.ylabel(f"Number of construct above {min_bases_cov} reads")
        plt.grid()

    def base_coverage_for_all_constructs(df, min_bases_cov):
        plt.figure(figsize=(20, 8))
        plt.plot(np.array(df['cov_bases_var'].sort_values(ascending=False).reset_index())[:,1])
        plt.plot(np.arange(0,int(df.construct.count()),1), [min_bases_cov]*int(df.construct.count()))
        plt.legend(['Dataframe', '1000 reads line'])
        plt.xlabel('Constructs (sorted)')
        plt.ylabel('# of reads of the worst covered base in the var part for a given structure in a given tube')

    def random_base_coverage_plot_wise(df, min_bases_cov):
        random_selection = np.random.randint(len(df), size=(9))
        fig = plt.figure(figsize=(25, 10))
        for i in range(9):
            axes1 = plt.subplot(int('33'+str(i+1)))
            plt.plot(np.array(df['cov_bases'].iloc[random_selection[i]]))
            start, end = df['var_start_index'].iloc[random_selection[i]], df['var_end_index'].iloc[random_selection[i]]
            plt.plot(np.arange(start, end, 1), np.array(df['cov_bases'].iloc[random_selection[i]])[start:end])
            plt.plot(np.arange(0, len(df['cov_bases'].iloc[random_selection[i]])), len(df['cov_bases'].iloc[random_selection[i]])*[min_bases_cov])
            plt.xlabel("Bases")
            plt.ylabel("Coverage")
            plt.title(f"Construct {df['construct'].iloc[random_selection[i]]}, tube {df['tube'].iloc[random_selection[i]]} ")
            plt.grid()
            plt.legend(["Base coverage (all)", 'Base coverage (var)', 'min_bases_cov'])
            axes2 = axes1.twinx()   # mirror them
            axes2.set_ylabel('Coverage [%]')
            axes2.set_ylim((0,100*max(df['cov_bases'].iloc[random_selection[i]]))/df['num_reads'].iloc[random_selection[i]])
        fig.tight_layout()

    def heatmap(df, column):
        base_cov_plot = df.pivot("tube","construct", column).astype(float)
        f, ax = plt.subplots(figsize=(28, 10))
        sns.heatmap(base_cov_plot, annot=False, linewidths=0, ax=ax, norm=LogNorm())


    def fit_deltaG(df, tube): #TODO
               # Fit
        fit = lambda a, b, c, dG: a/(1+b*np.exp(-dG/(R*T))) + c

            ## Max mut freq
        #a = np.array(df_use['mod_bases_A'].loc[tube][-1]+df_use['mod_bases_C'].loc[tube][-1]).mean()/np.array(df_use['info_bases'].loc[tube][-1]).mean()
            ## Baseline mut freq
        #   c = np.array(df_use['mod_bases_A'].loc[tube][0]+df_use['mod_bases_C'].loc[tube][0]).mean()/np.array(df_use['info_bases'].loc[tube][0]).mean()
            ## #TODO b
        # b = 1
        #  print(f"For tube {tube}, max mut freq a is {a}, b is {b}, baseline mut freq c is {c}")
        # plt.plot(t, fit(a, b, c, t))



    def mutation_rate(df, tube, construct, plot_type, index):
        
        df_use = df.set_index(['tube','construct'])
        
        if not plot_type in ['sequence','partition']:
            raise f"{plot_type} must be 'sequence' or 'partition', please check this argument"

        if plot_type == 'sequence':  # Plot the mutation rate for each base along the sequence

            mut_per_base = pd.DataFrame({'mut_rate': pd.Series(np.array(df_use[f"mut_bases"].loc[tube, construct][1:])/np.array(df_use[f"info_bases"].loc[tube, construct][1:]), dtype=object)
                                        ,'base':list(df_use['sequence'].loc[tube, construct])})\
                                        .reset_index()\
                                        .set_index(['base', 'index'])

            df_hist = pd.DataFrame()
            df_hist.index = mut_per_base.reset_index()['index']

            for base in ['A','C','G','T']:
                df_hist[base] = pd.Series(dtype=float)
                df_hist[base] = mut_per_base.loc[base]

            if index == 'base':
                df_hist.index = mut_per_base.reset_index()['base']

            ax = df_hist.plot.bar(stacked=True, figsize=(35,7))
            plt.title(f"tube {tube}, construct {construct}")

        if plot_type == 'partition': # Plot the partition of mutations for each base along the sequence
            df_hist = pd.DataFrame()
            for base in ['A','C','G','T']:
                df_hist[f"mod_bases_{base}"]  = np.array(df_use[f"mod_bases_{base}"].loc[tube, construct][1:])/df_use['info_bases'].loc[tube, construct][1:]

            if index == 'base':
                df_hist.index = list(df_use['sequence'].loc[tube,construct])

            ax = df_hist.plot.bar(stacked=True, figsize=(35,7))

    def deltaG(df, tube):
        df_use = df.set_index(['tube','construct'])

        fig = utils.define_figure(title=tube,
                                            xlabel='deltaG [kcal]',
                                            ylabel='Mutation ratio',
                                            figsize=(20,5))

        stack_for_plot = {'0':{'x':[],'y':[]},'1':{'x':[],'y':[]}}

        for construct in df.construct.unique():
            var_part = utils.get_var_info(df=df, tube=tube, construct=construct)
            for base in ['A','C']:
                for var_struct_comp in ['0','1']:
                    try:    
                        this_base_mut =  var_part.xs((base,True,var_struct_comp), level=('base','paired','var_structure_comparison'))
                        stack_for_plot[var_struct_comp]['x'].extend(this_base_mut['var_deltaG'].to_list())
                        stack_for_plot[var_struct_comp]['y'].extend(this_base_mut['mut_rate'].to_list())
                    except:
                        continue
        plt.plot(stack_for_plot['0']['x'],stack_for_plot['0']['y'],'b.')
        plt.plot(stack_for_plot['1']['x'],stack_for_plot['1']['y'],'r.')
        plot.fit_deltaG(df_use, tube)
        plt.legend(['A and C bases of the variable part, predicted paired by RNAstructure for both the var part sequence and the full sequence',\
                    'A and C bases of the variable part, predicted paired by RNAstructure for the full sequence but not for the variable part sequence'])
        plt.ylim([0,0.15])
        fig.tight_layout()

    def compare_2_tubes(df, tubes, construct, axs=None):
        if axs is None:
            fig, axs = plt.subplots(1,1)

        paired = {True: '.',False:'x'}
        var_structure_comparison = {'0':'b','1':'r'}
        x_all, y_all = [], []
        for is_paired in paired: 
            for vsc in var_structure_comparison:
                try:
                    x, y = np.array(utils.get_var_info(df, tubes[1], construct)['mut_rate'].xs((is_paired,vsc), level=('paired','var_structure_comparison')), dtype=float),\
                            np.array(utils.get_var_info(df, tubes[0], construct)['mut_rate'].xs((is_paired,vsc),level=('paired','var_structure_comparison')), dtype=float)
                    axs.plot(x,y,f"{var_structure_comparison[vsc]}{paired[is_paired]}")
                    x_all.extend(x), y_all.extend(y)
                except:
                    axs.plot()
                    continue
        p =  np.poly1d(np.polyfit(x_all,y_all,1))
        t = np.linspace(min(x_all),max(x_all))
        axs.plot(t,p(t),'g-')
        axs.set(xlabel=f"Mutation rate of tube {tubes[1]}", ylabel=f"Mutation rate of tube {tubes[0]}")


    def compare_n_tubes(df, tubes, construct):
        fig, axs = plt.subplots(len(tubes)-1,len(tubes)-1, figsize= (25,25), sharex=True, sharey=True)
        for x in range(0,len(tubes)):
            for y in range(x):
                if x != y:
                    plot.compare_2_tubes(df, (tubes[x], tubes[y]), construct, axs[x-1][y])
        axs[0,len(tubes)-2].plot(0,0,'b.',0,0,'r.',0,0,'bx',0,0,'rx',0,0,'g-')            
        axs[0,len(tubes)-2].legend(['Paired in full sequence RNAstructure, paired in var sequence RNAstructure',
                    'Paired in full sequence RNAstructure, not paired in var sequence RNAstructure',
                    'Not paired in full sequence RNAstructure, not paired in var sequence RNAstructure',
                    'Not paired in full sequence RNAstructure, paired in var sequence RNAstructure',
                    'Fit'])



    def save_fig(path,title):
        full_path = utils.make_path(path)
        plt.savefig(f"{full_path}/{title}")


class utils:
    def make_path(path):
        path = os.path.normpath(path)
        path=path.split(os.sep)
        path[path.index('date')] = str(datetime.datetime.now())[:10]
        full_path = ''
        for repo in path:
            full_path = full_path + f"{repo}/"
            if not exists(full_path):
                os.mkdir(full_path)
        return full_path

    def connect_to_firebase(verbose = True):
        with open('dreem-542b7-firebase-adminsdk-ejldl-510f9fa7da.json','w') as json_file:
            json.dump({
        "type": "service_account",
        "project_id": "dreem-542b7",
        "private_key_id": "510f9fa7da38e277c2e9673e6550426ae9183267",
        "private_key": "-----BEGIN PRIVATE KEY-----\nMIIEvAIBADANBgkqhkiG9w0BAQEFAASCBKYwggSiAgEAAoIBAQCgExvQjKxHwVuw\nt1O6t8R4IoAdpl7CqcvtxuMw/jIfGvAnITpzJoB+O+AeIxwlSBO1bRCmoCq7T8Uj\nkJdzaYxo9fT/UTFWc3Uo3M2fmGY5YgVQpzHSIdWDD9r2bwgwkXS9bPxk41yYADsH\n1b62lN62akrbNp6MTB2Vib+qWJqX/F34BhEDxGa7g31yPSk8fc7NyEgHHoii14Xl\ncV9JZ0KsqiDr9PJKJ6Jl8zpr9R9X3mCYi50vJYcRUKCJ6bxtWN933lHhnq7rusuo\nIdOUA83RTXM7RoGAQALQm5cyBPiO9xMW20OBVOGN/kfWsIs9/GQerRtLGAnRN0t0\ns9WWnpdHAgMBAAECggEADu6pZhNxUMJDTOFVILJa1AAX5mwqI8uWF+i5Mc1MnKU1\nKNlLLAm3686nEfihfALUv9RcPMbtJYsD71TiI+SBMhtbjuOikBd2IukyD0S2qHyx\n1Tu7hIgedDrq6Jkj8O/orXD4vGqPLSi8WPdB8qNBgU+6Cuf18014ZwYyCHB6f1no\nH3IMzsKQL10p/2vfA5PmYDOOhWekLV34sqiPerPYnuyiWmclgM4L+smpTOo1k1h8\neK1t/uf0SZQgQvXraibgTVs38MDd2mAo/JSicmVU+5PK22zdMfptvia2cfx2HGP/\n8ZRd8pDJi20lUgNr1yI9kUNDBPp/umVK6XwUftJjAQKBgQDLl7ckNupKV2d2aCdo\nKUs6NQZHiV1nRA6lbLQU51lzDKVhtTGLM75EojiUTPmOq9AWQGeHm9iUmv7UaLtg\nDeCB11AAbKFAxyxwjxE0E1d/ii8kh2ZYTAUh64viU/5dWLX+ZT2g8P9r92xOYEle\nt1bMs953qpFRXYf9GpCBCe0rAQKBgQDJR6hr1aarmuNAXLAgnSgn2CNEZHSUi7sT\nkvaG/EOl7XnLq4fMV5uwA1n/r/5pufn9KaS3MnJlUpqqcBtlmnAzvi5pvtKmpHC8\nQfimE15smUkgfCb8l+LesWk1lM/kpSeAn0L/FqL26KLnApa/4VKhLIFI+sPen/y4\nHZZ3FEmqRwKBgHUnSGu+bfN5eD/aj1KQ8Ij+Gi7wDJ9vuj3W34ln10Es9b3T1j6T\n99jmwEgWQ0Sl+YfUZ77RHz/kMN9ppOkREy+kBpU37VKpShk7OlsNBjyN97K9d1c3\n53wtXsFONADjG1bYSy5hf5lRNzGilpW6Smhg2JNjw1texvIOZzjZzXABAoGAYK0w\ncgr+sPIGMQXT+vZBMVIZLmJptGehBXfTPWaxP2Ne2rqa0UVLHDGf6rWnpzSSpEx6\nNxvd4ljYvQB3yEdzmQbB2Dy1hSD6nRG60ln/Qn4lp5q6RxzU9U2VUQ0XBaVl4dud\nHFTNFXcLt5WAvs0FGTD9MAZySd3iTrS3bp6p+0UCgYBnwQ3JTz5SyNkR7WqYlVgK\n2OIeslo7hCqGnn7pyuggNDr/YrE/P8i3cly+6BJ0twJMVi+qg1WqOKiy5AQhHt+5\nXk68ZQu1911T/saO76lafCGqZkdGlsi7k1uv3lh8/SyqXcXKGy9L1KJKTI7YXfbK\nvV2Ga5Uce124cBoekQc9Dw==\n-----END PRIVATE KEY-----\n",
        "client_email": "firebase-adminsdk-ejldl@dreem-542b7.iam.gserviceaccount.com",
        "client_id": "111801627848439891468",
        "auth_uri": "https://accounts.google.com/o/oauth2/auth",
        "token_uri": "https://oauth2.googleapis.com/token",
        "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
        "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/firebase-adminsdk-ejldl%40dreem-542b7.iam.gserviceaccount.com"
        },json_file )

        cred = credentials.Certificate("dreem-542b7-firebase-adminsdk-ejldl-510f9fa7da.json")
        # Initialize the app with a service account, granting admin privileges
        # As an admin, the app has access to read and write all data, regradless of Security Rules

        try:
            default_app = firebase_admin.initialize_app(cred, {
                'databaseURL':'https://dreem-542b7-default-rtdb.firebaseio.com/'
                })
        except:
            if verbose: print('Re-used the previous firebase connection')


    def define_figure(title, xlabel, ylabel, figsize):
        fig = plt.figure(figsize=figsize)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        return fig

    def added_content_per_construct(df, attribute):          
        return df.set_index('construct').sort_values(attribute)[attribute].groupby('construct').apply(lambda x:np.array(x)[0]).sort_values()


    def get_var_info(df, tube, construct):
        df_use = df.set_index(['tube','construct'])
        start, end = df_use['var_start_index'].loc[(tube,construct)] , df_use['var_end_index'].loc[(tube,construct)]     
        mut_per_base = pd.DataFrame({'mut_rate':pd.Series(np.array(df_use[f"mut_bases"].loc[tube, construct][1:])/np.array(df_use[f"info_bases"].loc[tube, construct][1:]), dtype=object),
                                'base':list(df_use['sequence'].loc[tube, construct]),
                                'paired': np.array([bool(x != '.') for x in list(df_use['full_structure'].loc[tube,construct])]),\
                                'var_structure_comparison': pd.Series(list(df_use['var_structure_comparison'].loc[tube,construct]), index=list(range(start, end)))\
                                ,'var_deltaG':df_use['var_deltaG'].loc[tube, construct]})\
                                .dropna()\
                                .reset_index()\
                                .set_index(['base', 'paired', 'var_structure_comparison','index'])
        return mut_per_base


    def columns_to_csv(df, tubes, columns, title, path):
        full_path = utils.make_path(path)
        df_print = df[df.tube.isin(tubes)]
        df_print = df_print[columns] 
        np.set_printoptions(suppress=True)
        df_print['mut_rate'] = df_print.apply(lambda row: np.float32(np.array(row['mut_bases'])/np.array(row['info_bases'])), axis=1)
        df_print.to_csv(f"{full_path}/{title}.csv")

    def deltaG_vs_construct_to_csv(df, title, path, tubes):
        full_path = utils.make_path(path)
        df[df['tube']==tubes[0]][['construct','var_deltaG','full_deltaG']].reset_index().drop(columns=['index']).to_csv(f"{full_path}/{title}")


    def big_script(df, tubes, constructs, analyse, study):
        show_plots = False
        if analyse['base_per_base_partition'] or analyse['base_per_base_sequence'] or analyse['deltaG']:
            print('base per base and/or deltaG processing')
        for tube in tubes:
            print(tube, end=' ')
            for construct in constructs:
                if analyse['base_per_base_partition']:
                    plot.mutation_rate(df, tube, construct, plot_type='partition', index='base')  
                    plot.save_fig(path=f"data/figs/date/{study}/mut_per_base/partition/{construct}/", 
                                title=f"base_per_base_partition_{tube}_{construct}")
                    plt.close(not show_plots)
                if analyse['base_per_base_sequence']:
                    plot.mutation_rate(df, tube, construct, plot_type='sequence', index='index')
                    plot.save_fig(path=f"data/figs/date/{study}/mut_per_base/sequence/{construct}/", 
                                title=f"base_per_base_partition_{tube}_{construct}")
                    plt.close(not show_plots)
            if analyse['deltaG']:
                plot.deltaG(df=df, tube=tube)
                plot.save_fig(path=f"data/figs/date/{study}/deltaG/", 
                    title=f"deltaG_{tube}")
                plt.close(not show_plots)

        if analyse['tube_comparison']:
            print('Tube comparison processing')
            for construct in constructs:
                    plot.compare_n_tubes(df, tubes, construct)
                    plot.save_fig(path=f"data/figs/date/comparison/{study}", 
                                title=f"comparison_{study}_{construct}")
                    plt.close(not show_plots)
                    print(construct, end=' ')
        if analyse['columns_csv']:
            print('Push columns to csv')
            utils.deltaG_vs_construct_to_csv(df=df, title=f"deltaG_vs_construct", path = f"data/figs/date", tubes=tubes)

        if analyse['deltaG_construct']:
            print('Push deltaG vs construct to csv')
            utils.columns_to_csv(df=df,
                            tubes=tubes,
                            columns=['tube', 'construct','full_sequence','var_sequence','mut_bases','info_bases'],
                            title=f"about_{study}",
                            path=f"data/figs/date/{study}"
                            )

    def deltaG_vs_construct_to_csv(df, title, path, tubes):
        full_path = utils.make_path(path)
        df[df['tube']==tubes[0]][['construct','var_deltaG','full_deltaG']].reset_index().drop(columns=['index']).to_csv(f"{full_path}/{title}")
        