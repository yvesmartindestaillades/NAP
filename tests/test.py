from dreem_nap import utils
from dreem_nap.study import Study

def test_figs(temp:Study):
    pl = ["mut_histogram('A1','9572','index')"
    #'mut_rate_vs_base_non_pairing_prob()',
    #'sliding_window_r2_gini()',
   # 'study_base()',
   # 'study_sample()',
   # 'base_wise_mut_vs_prob()',
   # 'deltaG()',
   # 'deltaG_basewise()',
  #  'correlation_n_samples()',
  #  'base_coverage()',
   # 'base_coverage_ROI_for_all_constructs()',
  #  'random_9_base_coverage()',
  #  'sample_coverage_distribution()',
 #   'valid_construct_per_sample()',
   # 'heatmap()'
   ]

    for p in pl:
        eval(f"temp.{p}()")
        utils.save_fig(f"test/figs/{pl[:-2]}.png")
        print(f"passed {p}")

if __name__ == '__main__':
    temp = Study('temperature',['A1','B2','B3'], [10, 20, 30], 'Example values [no unit]', 'Just an example study')
    temp.load_df_from_local_files('data/DEMULTIPLEXED',10)
    test_figs(temp)