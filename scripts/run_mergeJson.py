#!/usr/bin/python
import os

version = 'v1'
directory="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonNtuple_json/Run3/V1p19/Data2022/{0}/".format(version)
golden_json_path = '../data/Certification/20221215/'
samples = [
'DisplacedJet-EXOCSCCluster_Run2022E-PromptReco-{0}'.format(version),
'DisplacedJet-EXOCSCCluster_Run2022F-PromptReco-{0}'.format(version),
'DisplacedJet-EXOCSCCluster_Run2022G-PromptReco-{0}'.format(version),
]
if not os.path.exists(directory+'merged'): 
	print(directory + 'merged')
	os.makedirs(directory + 'merged')
for sample in samples:
	print(sample)
	input_json = '{}/{}/*.json'.format(directory, sample)
	output_json = '{}/merged/{}.json'.format(directory, sample)
	goodLumi_json = '{}/merged/{}_goodLumi.json'.format(directory, sample)
	if '2022E' in sample:cert = golden_json_path + 'Cert_Collisions2022_eraE_359022_360331_Golden.json'
	elif '2022F' in sample: cert = golden_json_path + 'Cert_Collisions2022_eraF_360390_362167_Golden.json'
	elif '2022G' in sample: cert = golden_json_path + 'Cert_Collisions2022_eraG_362433_362760_Golden.json' 
	else: assert(False)
	merge_cmd = 'mergeJSON.py {} --output {}'.format(input_json, output_json)
	filter_cmd = 'compareJSON.py --and {}  {} > {}'.format(output_json, cert, goodLumi_json)
	print(merge_cmd)
	os.system(merge_cmd)
	os.system(filter_cmd)
