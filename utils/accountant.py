    # Data structure as a JSON
    # 	[Tube 1 name]:
    #       (log content, example here)
	#		Sample: 210422Rou_D21-3900
	#		Reference_file: /lab/solexa_rouskin/projects/mfallan/rna_np_210422Rou/Ref_Genome/23S.fasta
	#		Reference_genomes: ['23S']
	#		Paired: True
	#		Num_surrounding bases_for_del: 10
	#		Q_score_file: /lab/solexa_rouskin/phred_33.txt
	#		Q_score_cutoff: 20
	#		Time_taken: 0.63 mins
	#		Finished_at: 2021-05-09 00:10

	#		(quality control content, example here)
	#		Signal_to_noise_ratio: 0.9
	#		#TODO

	#		DMS_mutation:
	# 			Position:	
	# 				Mismatches	
	# 				Mismatches + Deletions
	#			...

	#		Read_coverage:
	#			Position:
	#				Read_coverage_value
	#	[Tube 2 name]
	#		...
	#   ...

import json, os

from matplotlib.font_manager import json_load

JSON_LOCATION = "data/JSON output/"
BITVECTOR_LOCATION = "data/BitVector_Plots/"

class Accountant():
	def __init__(self, ref_gen_coord_wrt_tubes= None):
		# Dict of the reference genome w.r.t samples (tubes) 
		self.ref_gen_coord_wrt_tubes = ref_gen_coord_wrt_tubes 
	
	def fetch_ref_gen_coord_wrt_tubes():
		# Fetches reference genomes and samples names
		# to read by reading the files names
		return None

	def log_content_to_dict(self, file_prefix):
		# Turns a log file into a dictionary
		with open(BITVECTOR_LOCATION + file_prefix+'_log.txt') as f:
			log_content = {}
			lines = f.readlines()
			attributes = ["Sample", "Reference file","Reference genomes","Paired"\
						 ,"Num surrounding bases for del","Q score file","Q score cutoff",\
						 "Time taken","Finished at"]
			for attr in attributes:
				for line in lines:
					if attr in line:
						start = line.index(':')+2
						log_content[attr] = line[start:]
		return None
	
	def popavg_react_txt_to_dict(self, file_prefix):
		# Turns a txt popavg_react file into a dictionary
		return None

	def read_coverage_html_to_dict(self, file_prefix):
		# Turns a read coverage html file into a dictionary
		return None

	def make_JSON(self, file_name):
		# Reads every sample one-by-one and dumps it into a json file
		with open(JSON_LOCATION+file_name, 'w') as file_object:  #open the file in write mode
			for tube in self.ref_gen_coord_wrt_tubes.keys():
				data_tube = {}
				data_tube[tube] = self.log_content_to_dict(file_prefix=tube)
				file_prefix = tube+'_'+self.ref_gen_coord_wrt_tubes[tube]+'_'
				data_tube[tube]['DMS_mutation'] = self.popavg_react_txt_to_dict(file_prefix=file_prefix)
				data_tube[tube]['Read_coverage'] = self.read_coverage_html_to_dict(file_prefix=file_prefix)
				json.dump(data_tube, file_object)


# Validation script
if __name__ == "__main__":
	ref_gen_coord_wrt_tubes = {
		'210422Rou_D21-3897_23S':'25_463',
		'210422Rou_D21-3900_rsc1218v1n924':'436_892'
	}
	accountant = Accountant(ref_gen_coord_wrt_tubes)
	accountant.make_JSON('my_JSON.json')