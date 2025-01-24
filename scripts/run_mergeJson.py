#!/usr/bin/python
import os

version = 'v4'
directory="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonNtuple_json/Run3/V1p19/2024/{0}/".format(version)
golden_json_path = '../data/Certification/'
if not os.path.exists(directory+'merged'): 
	print(directory + 'merged')
	os.makedirs(directory + 'merged')
samples = os.listdir(directory)
for sample in samples:
	print(sample)
	if "merge" in sample:continue
	input_json = '{}/{}/*.json'.format(directory, sample)
	output_json = '{}/merged/{}.json'.format(directory, sample)
	goodLumi_json = '{}/merged/{}_goodLumi.json'.format(directory, sample)
	if '2022' in sample:cert = golden_json_path + 'Cert_Collisions2022_355100_362760_Golden.json'
	elif '2023' in sample: cert = golden_json_path + 'Cert_Collisions2023_366442_370790_Golden.json'
	elif '2024' in sample: cert = golden_json_path + 'Cert_Collisions2024_378981_386951_Golden.json'
	else: assert(False)
	merge_cmd = 'mergeJSON.py {} --output {}'.format(input_json, output_json)
	filter_cmd = 'compareJSON.py --and {}  {} > {}'.format(output_json, cert, goodLumi_json)
	print(merge_cmd)
	os.system(merge_cmd)
	os.system(filter_cmd)

