import os
#mc
version = '/V1p19/'
base_path  = "/eos/uscms//store/group/lpclonglived/displacedJetMuonNtuple/" + version 
base_list_path = "../lists/displacedJetMuonNtuple/" + version
year = ['MC_Summer22EE', 'MC_Summer22','MC_Summer23BPix', 'MC_Summer23']


for y in year:
    path = base_path + y + '/v2/'
    list_path = base_list_path + y + '/v2/'
    os.system("mkdir -p " + list_path)
    samples = os.listdir(path)
    for s in samples:
        if not "DY" in s:continue
        with open(list_path + s + '.txt', 'w') as fp:
                    path_temp = path + s + '/'
                    for root, dirs, files in os.walk(path_temp):
                            for filename in files:
                                    fp.write("root://cmsxrootd.fnal.gov//" + os.path.join(root, filename).replace("/eos/uscms/","") + "\n")
