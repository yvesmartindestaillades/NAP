=========================
Dataframe columns content
=========================

Here, you'll learn how to build a database that's compatible with NAP. 

Each row of the csv file is a **mutation profile**, which is a unique
combination of the columns “sample”, “construct”, “section” and
“cluster”. Each column of the csv is an attrbute of the mutation
profile. Attributes comes from different files, such as ``library.csv``,
``samples.csv``, and by running softwares such as ``DREEM``,
``RNAstructure`` and ``Poisson``.

samples.csv
========================

+-----------------+-----------------+-----------------+-----------------+
| attribute       | type            | description     | comment         |
+=================+=================+=================+=================+
| sample          | str             | fastq file      | cannot be int   |
|                 |                 | prefix          |                 |
+-----------------+-----------------+-----------------+-----------------+
| user            | str             | Who did the     |                 |
|                 |                 | experiment      |                 |
+-----------------+-----------------+-----------------+-----------------+
| date            | str             | Date of the     |                 |
|                 |                 | experiment      |                 |
+-----------------+-----------------+-----------------+-----------------+
| exp_env         | str             | Experimental    | Can only be one |
|                 |                 | environment,    | of the two      |
|                 |                 | “in_vivo” or    | options         |
|                 |                 | “in_vitro”      | “in_vivo” and   |
|                 |                 |                 | “in_vitro”      |
+-----------------+-----------------+-----------------+-----------------+
| temperature_k   | float           | Temperature in  |                 |
|                 |                 | Kelvin          |                 |
+-----------------+-----------------+-----------------+-----------------+
| in              | float           | Total           |                 |
| c_time_tot_secs |                 | incubation time |                 |
|                 |                 | in seconds      |                 |
+-----------------+-----------------+-----------------+-----------------+
| buffer          | str             | Exact buffer    | Only if exp_env |
|                 |                 | including Mg,   | == “in_vitro”   |
|                 |                 | eg 300mM Sodium |                 |
|                 |                 | Cacodylate ,    |                 |
|                 |                 | 3mM Mg          |                 |
+-----------------+-----------------+-----------------+-----------------+
| cell_line       | str             | Cell line       | Only if exp_env |
|                 |                 |                 | == “in_vivo”    |
+-----------------+-----------------+-----------------+-----------------+

library.csv
========================

+-----------------+-----------------+-----------------+-----------------+
| attribute       | type            | description     | comment         |
+=================+=================+=================+=================+
| construct       | str             | fasta file      | cannot be int   |
|                 |                 | constructs      |                 |
|                 |                 | names           |                 |
+-----------------+-----------------+-----------------+-----------------+
| section         | str(list(       | sub-sequences   | include by      |
|                 | {str:tuple(int, | of a construct  | default         |
|                 | int)}))         |                 | {‘ful           |
|                 |                 |                 | l’:(min(index), |
|                 |                 |                 | max(index))}    |
+-----------------+-----------------+-----------------+-----------------+
| [attribute]     | [type]          | a per-construct | can be family   |
|                 |                 | attribute       | or flank for    |
|                 |                 |                 | example         |
+-----------------+-----------------+-----------------+-----------------+

RNAstructure
========================

+-----------------+-----------------+-----------------+-----------------+
| attribute       | type            | description     | comment         |
+=================+=================+=================+=================+
| deltaG_min      | float           | minimum energy  |                 |
|                 |                 | for the         |                 |
|                 |                 | sequence        |                 |
+-----------------+-----------------+-----------------+-----------------+
| deltaG_min_T    | float           | minimum energy  |                 |
|                 |                 | for the         |                 |
|                 |                 | sequence using  |                 |
|                 |                 | temperature     |                 |
+-----------------+-----------------+-----------------+-----------------+
| deltaG_min_DMS  | float           | minimum energy  |                 |
|                 |                 | for the         |                 |
|                 |                 | sequence and    |                 |
|                 |                 | DMS signal      |                 |
+-----------------+-----------------+-----------------+-----------------+
| d               | float           | minimum energy  |                 |
| eltaG_min_T_DMS |                 | for the         |                 |
|                 |                 | sequence using  |                 |
|                 |                 | temperature and |                 |
|                 |                 | DMS signal      |                 |
+-----------------+-----------------+-----------------+-----------------+
| deltaG_ens      | float           | average energy  |                 |
|                 |                 | of the          |                 |
|                 |                 | partition       |                 |
|                 |                 | function for    |                 |
|                 |                 | this sequence   |                 |
+-----------------+-----------------+-----------------+-----------------+
| structure       | str             | minimum energy  |                 |
|                 |                 | structure for   |                 |
|                 |                 | the sequence    |                 |
+-----------------+-----------------+-----------------+-----------------+
| structure_T     | str             | structure for   |                 |
|                 |                 | the sequence    |                 |
|                 |                 | using           |                 |
|                 |                 | temperature as  |                 |
|                 |                 | an input        |                 |
+-----------------+-----------------+-----------------+-----------------+
| structure_DMS   | str             | energy          |                 |
|                 |                 | structure for   |                 |
|                 |                 | the sequence    |                 |
|                 |                 | using DMS as an |                 |
|                 |                 | input           |                 |
+-----------------+-----------------+-----------------+-----------------+
| structure_T_DMS | str             | energy          |                 |
|                 |                 | structure for   |                 |
|                 |                 | the sequence    |                 |
|                 |                 | using           |                 |
|                 |                 | temperature and |                 |
|                 |                 | DMS as an input |                 |
+-----------------+-----------------+-----------------+-----------------+

Poisson
========================

+-----------------+-----------------+-----------------+-----------------+
| attribute       | type            | description     | comment         |
+=================+=================+=================+=================+
| poisson_min     | s               | Low boundary of |                 |
|                 | tr(list(float)) | the Poisson     |                 |
|                 |                 | confidence      |                 |
|                 |                 | interval for    |                 |
|                 |                 | the mutation    |                 |
|                 |                 | rate            |                 |
+-----------------+-----------------+-----------------+-----------------+
| poisson_max     | s               | High boundary   |                 |
|                 | tr(list(float)) | of the Poisson  |                 |
|                 |                 | confidence      |                 |
|                 |                 | interval for    |                 |
|                 |                 | the mutation    |                 |
|                 |                 | rate            |                 |
+-----------------+-----------------+-----------------+-----------------+
| poisson_low     | s               | Length of the   | mut_r           |
|                 | tr(list(float)) | low error bar   | ate-poisson_min |
|                 |                 | for mutation    |                 |
|                 |                 | rate            |                 |
+-----------------+-----------------+-----------------+-----------------+
| poisson_high    | s               | Length of the   | poiss           |
|                 | tr(list(float)) | high error bar  | on_max-mut_rate |
|                 |                 | for mutation    |                 |
|                 |                 | rate            |                 |
+-----------------+-----------------+-----------------+-----------------+

DREEM
========================

+-----------------+-----------------+-----------------+-----------------+
| attribute       | type            | description     | comment         |
+=================+=================+=================+=================+
| sample          | str             | fastq file      | cannot be int   |
|                 |                 | prefix          |                 |
+-----------------+-----------------+-----------------+-----------------+
| construct       | str             | fasta file      | cannot be int   |
|                 |                 | constructs      |                 |
|                 |                 | names           |                 |
+-----------------+-----------------+-----------------+-----------------+
| cluster         | str             | alternative     | default:        |
|                 |                 | mutational      | ‘pop_avg’       |
|                 |                 | profiles given  |                 |
|                 |                 | by DREEM        |                 |
+-----------------+-----------------+-----------------+-----------------+
| sequence        | str             | nucleotides     | uses A, C, G, T |
|                 |                 | sequence        |                 |
+-----------------+-----------------+-----------------+-----------------+
| data_type       | str             | Type of data    | useful?         |
|                 |                 | (DMS)           |                 |
+-----------------+-----------------+-----------------+-----------------+
| num_reads       | int             | Number of reads | useful?         |
|                 |                 | for this        |                 |
|                 |                 | mutation        |                 |
|                 |                 | profile         |                 |
+-----------------+-----------------+-----------------+-----------------+
| num_aligned     | int             |                 | useful?         |
+-----------------+-----------------+-----------------+-----------------+
| n               | str(list(int))  | Count of        | useful?         |
| um_of_mutations |                 | mutations per   |                 |
|                 |                 | read            |                 |
+-----------------+-----------------+-----------------+-----------------+
| mut_bases       | str(list(int))  | Per-residue     | 0-indexed       |
|                 |                 | count of        |                 |
|                 |                 | mutations       |                 |
+-----------------+-----------------+-----------------+-----------------+
| info_bases      | str(list(int))  | Per-residue     | Diff with cov   |
|                 |                 | count of valid  | bases?          |
|                 |                 | reads           |                 |
+-----------------+-----------------+-----------------+-----------------+
| cov_bases       | str(list(int))  | Per-residue     |                 |
|                 |                 | count of        |                 |
|                 |                 | covered bases   |                 |
+-----------------+-----------------+-----------------+-----------------+
| del_bases       | str(list(int))  | Per-residue     | useful?         |
|                 |                 | count of        |                 |
|                 |                 | deleted bases   |                 |
+-----------------+-----------------+-----------------+-----------------+
| ins_bases       | str(list(int))  | Per-residue     | useful?         |
|                 |                 | count of        |                 |
|                 |                 | inserted bases  |                 |
+-----------------+-----------------+-----------------+-----------------+
| mod_bases_A     | str(list(int))  | Per-residue     | useful?         |
|                 |                 | count of        |                 |
|                 |                 | mutations to a  |                 |
|                 |                 | A base          |                 |
+-----------------+-----------------+-----------------+-----------------+
| mod_bases_C     | str(list(int))  | Per-residue     | useful?         |
|                 |                 | count of        |                 |
|                 |                 | mutations to a  |                 |
|                 |                 | C base          |                 |
+-----------------+-----------------+-----------------+-----------------+
| mod_bases_G     | str(list(int))  | Per-residue     | useful?         |
|                 |                 | count of        |                 |
|                 |                 | mutations to a  |                 |
|                 |                 | G base          |                 |
+-----------------+-----------------+-----------------+-----------------+
| mod_bases_T     | str(list(int))  | Per-residue     | useful?         |
|                 |                 | count of        |                 |
|                 |                 | mutations to a  |                 |
|                 |                 | T base          |                 |
+-----------------+-----------------+-----------------+-----------------+
| mut_rates       | s               | Per-residue     | mut_ba          |
|                 | tr(list(float)) | count of        | ses/info_bases, |
|                 |                 | mutation        | shall we use    |
|                 |                 | divided by the  | cov_bases       |
|                 |                 | count of valid  | instead?        |
|                 |                 | reads           |                 |
+-----------------+-----------------+-----------------+-----------------+
| worst_cov_bases | int             | min(info_bases) | to adapt to     |
|                 |                 | (or cov_bases?) | per-section and |
|                 |                 |                 | per-cluster     |
|                 |                 |                 | samples         |
+-----------------+-----------------+-----------------+-----------------+
| s               | int             | number of reads | useful?         |
| kips_short_read |                 | that we don’t   |                 |
|                 |                 | use because     |                 |
|                 |                 | they are too    |                 |
|                 |                 | short.          |                 |
+-----------------+-----------------+-----------------+-----------------+
| skip            | int             | number of reads | useful?         |
| s_too_many_muts |                 | that that we    |                 |
|                 |                 | don’t use       |                 |
|                 |                 | because they    |                 |
|                 |                 | have so many    |                 |
|                 |                 | mutations, and  |                 |
|                 |                 | therefore we    |                 |
|                 |                 | have low        |                 |
|                 |                 | confidence.     |                 |
+-----------------+-----------------+-----------------+-----------------+
| skips_low_mapq  | int             | number of reads | useful?         |
|                 |                 | that that we    |                 |
|                 |                 | don’t use       |                 |
|                 |                 | because the map |                 |
|                 |                 | score is too    |                 |
|                 |                 | low (default is |                 |
|                 |                 | below 15)       |                 |
+-----------------+-----------------+-----------------+-----------------+
|                 |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
