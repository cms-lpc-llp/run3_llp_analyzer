import os
#mc
version = '/V1p19/'
base_path  = "/eos/uscms//store/group/lpclonglived/displacedJetMuonNtuple/" + version 
base_list_path = "../lists/displacedJetMuonNtuple/" + version
year = ['Data2024']
dataset_list=['JetMET0','JetMET1']
for y in year:
    for dataset in dataset_list:
        path = base_path + y + '/v1/' + dataset +'/'
        list_path = base_list_path + y + '/v1/'+ dataset +'/'
        os.system("mkdir -p " + list_path)
        samples = os.listdir(path)
        for s in samples:
            print(list_path + s + '.txt')
            with open(list_path + s + '.txt', 'w') as fp:
                        path_temp = path + s + '/'
                        for root, dirs, files in os.walk(path_temp):
                                for filename in files:
                                        fp.write("root://cmsxrootd.fnal.gov//" + os.path.join(root, filename).replace("/eos/uscms/","") + "\n")
