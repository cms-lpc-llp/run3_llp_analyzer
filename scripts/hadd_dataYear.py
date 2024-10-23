#!/bin/env python
import os
import glob
base_path = '/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/'
version = '/v16/normalized/'
output_path = base_path + 'Data_all' + version
os.system('mkdir -p {}'.format(output_path))
print(output_path)
datasetList = {}
file_name = {
        '2022': 'DisplacedJet-EXOCSCCluster_Run2022-PromptReco_goodLumi.root',
        '2023': 'Muon-EXOCSCCluster_Run2023-PromptReco_goodLumi.root',
        '2024': 'Muon-Run2024-PromptReco_goodLumi.root',
}
year = list(file_name.keys())

os.system("eval `scram runtime -sh`")
os.system("rm -f {}".format(file_name))

#file_name=HV_params_${portal}_m_${m}_ctau_${ctau}mm_${scale}.root
outputRoot= 'EXOCSCCluster_Run2022_2024_goodLumi.root'
cmd = f'hadd {outputRoot} '
cmd += ' '.join([f'{base_path}/Data{y}/{version}/{name}' for y,name in file_name.items()])
os.system(cmd)
if os.path.isfile(outputRoot):
        os.system('mv {0} {1}/{0}'.format(outputRoot, output_path))
        if os.path.isfile(output_path + outputRoot): print('copy succeeded')

