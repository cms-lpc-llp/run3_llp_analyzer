import os
#MC
path  = "/eos/uscms//store/group/lpchscp/ntuples/V1p0/MC_UL18/"
list_path = "lists/ntuples/V1p0/MC_UL18/"

#Data
path  = "/eos/uscms//store/group/lpchscp/ntuples/V1p0/Data_UL/SingleMuon/"
list_path = "lists/ntuples/V1p0/Data_UL/SingleMuon/"

os.system("mkdir -p " + list_path)
samples = os.listdir(path)
for s in samples:
        with open(list_path + s + '.txt', 'w') as fp:
                path_temp = path + s + '/'
                for root, dirs, files in os.walk(path_temp):
                        for filename in files:
                                fp.write("root://cmsxrootd.fnal.gov//" + os.path.join(root, filename).replace("/eos/uscms/","") + "\n")
