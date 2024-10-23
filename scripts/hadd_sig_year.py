#!/bin/env python
import os
import glob
base_path = '/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/'
version = '/v16/normalized/'
output_path = base_path + 'MC_all' + version
os.system('mkdir -p {}'.format(output_path))
print(output_path)
datasetList = {}
year = ['Summer22EE','Summer23','Summer23BPix']
input_path = base_path + '/MC_Summer22EE/' + version
files = os.listdir(input_path)
for f in files:
        file_name  = '_'.join(f.split('_')[:-2])
        print(file_name)
        os.system("eval `scram runtime -sh`")
        os.system("rm -f {}".format(file_name))

        #file_name=HV_params_${portal}_m_${m}_ctau_${ctau}mm_${scale}.root
        outputRoot=file_name + '_50000pb_weighted.root'
        print(file_name, outputRoot)
        cmd = f'hadd {outputRoot} '
        cmd += ' '.join([f'{base_path}/MC_{y}/{version}/{file_name}*' for y in year])
        os.system(cmd)
        if os.path.isfile(outputRoot):
                os.system('mv {0} {1}/{0}'.format(outputRoot, output_path))
                if os.path.isfile(output_path + outputRoot): print('copy succeeded')

