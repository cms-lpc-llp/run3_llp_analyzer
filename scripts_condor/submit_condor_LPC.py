#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(
    description="Create and submit (or dry-run) a Condor workflow for run3_llp_analyzer on LPC."
)
parser.add_argument("sample", help="Dataset/sample name (must have a corresponding <sample>.txt in the list directory).")
parser.add_argument(
    "analyzer_str",
    help="Analyzer selector string (e.g. TnP_mdsnano).",
    choices=[
        "TnP",
        "TnP_mdsnano",
        "TnP_noClusters",
        "TrigEff",
        "TrigEff_mdsnano",
        "TnP_noClusters_VetoEff",
        "TnP_noClusters_VetoEff_mdsnano",
        "mdsnano",
        "mdsnano_hnl",
    ],
)
parser.add_argument(
    "output_end",
    help="Suffix appended to /store/group/lpclonglived/amalbert/HNL_Tau_Search/<output_end>/ for output.",
)
parser.add_argument("log_name", help="Tag used for local condor submit/log directories.")
parser.add_argument(
    "--dryRun",
    action="store_true",
    help="Do everything except eosmkdir + condor_submit (still writes the JDL).",
)
args = parser.parse_args()

sample = args.sample  # e.g. Muon0-Run2024B-PromptReco-v1-AOD (or a full path to a .txt list)
analyzer_str = args.analyzer_str  # e.g. TnP_mdsnano
output_end = args.output_end  # e.g. 2024_Data
log_name = args.log_name  # e.g. 2024_Data_Test

log_dir = "condor_job_output/log_{}".format(log_name)
submit_dir = "condor_job_output/submit_{}".format(log_name)
os.system("mkdir -p {}".format(log_dir))
os.system("mkdir -p {}".format(submit_dir))
executable = "scripts_condor/runAnalyzer_LPC.sh"
#analyzer = 'llp_MuonSystem_CA_TnP'
filesPerJob = 10
ntupler_version = 'V1p19/Data2023/'

#ntupler_version = "V1p19/MC_Summer22EE/v1/sixie/"


if analyzer_str=="TnP":analyzer="llp_MuonSystem_CA_TnP"
elif analyzer_str=="TnP_mdsnano":analyzer="llp_MuonSystem_CA_TnP_mdsnano"
elif analyzer_str=="TnP_noClusters":analyzer="llp_MuonSystem_CA_TnP_noClusters"
elif analyzer_str=="TrigEff":analyzer="llp_MuonSystem_CA_TrigEff"
elif analyzer_str=="TrigEff_mdsnano":analyzer="llp_MuonSystem_CA_TrigEff_mdsnano"
elif analyzer_str=="TnP_noClusters_VetoEff":analyzer="llp_MuonSystem_CA_TnP_noClusters_VetoEff"
elif analyzer_str=="TnP_noClusters_VetoEff_mdsnano":analyzer="llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano"
elif analyzer_str=="mdsnano":analyzer="llp_MuonSystem_CA_mdsnano"
elif analyzer_str=="mdsnano_hnl":analyzer="llp_MuonSystem_CA_mdsnano_hnl"
else:print("analyzer string not recognized, exiting ...");exit()
analyzer_version = 'v11'
#outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/{0}/{1}/".format(ntupler_version, analyzer_version)
#outputDirectoryBase="/store/group/lpclonglived/amalbert/Data_MC_Comp_TnP/results_from_cache_noSkim/{}/".format(output_end)
EOS_USER = os.getenv("USER") or "tlee"
outputDirectoryBase="/store/user/{}/HNL_Tau_Search/{}/".format(EOS_USER, output_end)
HOME = os.getenv("HOME") or ""
CMSSW_BASE = os.getenv("CMSSW_BASE")

# Prefer CMSSW_BASE if available, otherwise fall back to the repo location.
if CMSSW_BASE:
    Analyzer_DIR = os.path.join(CMSSW_BASE, "src", "run3_llp_analyzer")
else:
    Analyzer_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    print(f"[warn] CMSSW_BASE is not set; using Analyzer_DIR={Analyzer_DIR}")

if not Analyzer_DIR.endswith("/"):
    Analyzer_DIR += "/"

lists_root = os.path.join(Analyzer_DIR, "lists", "MDSNano", "v2")

def resolve_inputfilelist(sample_arg: str) -> str:
    """
    Allow either:
      - full path to a .txt list file
      - sample basename (without .txt), which will be searched under lists_root
    """
    if sample_arg.endswith(".txt"):
        candidate = sample_arg
        if not os.path.isabs(candidate):
            candidate = os.path.abspath(os.path.join(os.getcwd(), candidate))
        if os.path.exists(candidate):
            return candidate
        raise FileNotFoundError(candidate)

    matches = glob.glob(os.path.join(lists_root, "**", sample_arg + ".txt"), recursive=True)
    if len(matches) == 1:
        return matches[0]
    if len(matches) == 0:
        raise FileNotFoundError(f"No list file found for sample '{sample_arg}' under {lists_root}")
    raise RuntimeError(
        "Multiple list files matched sample '{0}'. Please pass the full .txt path.\n{1}".format(
            sample_arg, "\n".join(matches)
        )
    )

inputfilelist = resolve_inputfilelist(sample)
sample_name = os.path.basename(inputfilelist).replace(".txt", "")

try:
    rel = os.path.relpath(inputfilelist, lists_root)
    list_top = rel.split(os.sep, 1)[0]
except Exception:
    list_top = ""

is_data = list_top.startswith("Data") or "Run" in sample_name
isData = "--isData" if is_data else "no"


'''
#### Run on all signal samples ####
## year/isData
datasetList = OrderedDict()
samples = os.listdir(datasetListDir)
for s in samples:
    if not "EXOCSCCluster" in s:continue 
    if 'Data' in ntupler_version: datasetList[s.replace('.txt', '')] = ["2018", "yes"]
    else: datasetList[s.replace('.txt', '')] = ["2018", "no"]
############
'''



#for sample in datasetList.keys():

print("Preparing analyzer workflow for dataset: " + sample_name + "\n")

if not os.path.exists(inputfilelist):
    print("listfile: " + inputfilelist + " does not exist. exiting.")
    sys.exit(2)
cmd = ["awk", "END{print NR}",inputfilelist]
nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
maxjob=int(nfiles/filesPerJob)+1
mod=int(nfiles%filesPerJob)
if mod == 0: maxjob -= 1
outputDirectory = outputDirectoryBase + sample_name + "/"
if not args.dryRun:
    os.system(f"eos root://cmseos.fnal.gov mkdir -p {outputDirectory}")
else:
    print(f"[dryRun] Would run: eosmkdir {outputDirectory}")
print(sample_name)
if "Run2022E" in sample or "Run2022F" in sample or "Run2022G" in sample or "MC_Summer22EE" in sample:
    analyzerTag = "Summer22EE"
    jetVetoMap = "Summer22EE_23Sep2023_RunEFG_v1.root"
    pileupWeights = "PileupReweight_Summer22EE.root"
    HMT = "L1_efficiencies_2022_2023_032625-Hists-TEff.root"
    MET = "METTriggerEff_Summer22.root"
elif "Run2022C" in sample or "Run2022D" in sample or "MC_Summer22" in sample:
    analyzerTag = "Summer22"
    jetVetoMap = "Summer22_23Sep2023_RunCD_v1.root"
    pileupWeights = "PileupReweight_Summer22.root"
    HMT = "L1_efficiencies_2022_2023_032625-Hists-TEff.root"
    MET = "METTriggerEff_Summer22EE.root"
elif "Run2023D" in sample or "MC_Summer23BPix" in sample:
    analyzerTag = "Summer23BPix"
    jetVetoMap = "Summer23BPixPrompt23_RunD_v1.root"
    pileupWeights = "PileupReweight_Summer23BPix.root"
    HMT = "L1_efficiencies_2022_2023_032625-Hists-TEff.root"
    MET = "METTriggerEff_Summer23.root"
elif "Run2023B" in sample or "Run2023C" in sample or "MC_Summer23" in sample:
    analyzerTag = "Summer23"
    jetVetoMap = "Summer23Prompt23_RunC_v1.root"
    pileupWeights = "PileupReweight_Summer23.root"
    HMT = "L1_efficiencies_2022_2023_032625-Hists-TEff.root"
    MET = "METTriggerEff_Summer23BPix.root"
elif "Run2024" in sample_name or "MC_Summer24" in sample_name or "HNL"in sample_name or list_top in ("MC_Summer24", "Data2024"):
    analyzerTag = "Summer24"
    jetVetoMap = "Winter24Prompt24_2024BCDEFGHI.root"
    pileupWeights = "PileupReweight_Summer24.root"
    HMT = "HMT_Efficiencies_2024.root"
    MET = "METTriggerEff_Summer24.root"
else:
    print("Couldn't find a valid analyzer tag. Exiting ...")
    exit()


#####################################
#Create Condor JDL file
#####################################
#os.system("rm -f submit/{}_{}_Job*.jdl".format(analyzer, sample))
#os.system("rm -f log/{}_{}_Job*".format(analyzer, sample))

#Need to transfer in OPENSSL Libs and Bin for xrdcp to work in condor


jdl_file="{}/{}_{}_{}.jdl".format(submit_dir, analyzer, sample_name.replace("/","_"), maxjob)

tmpCondorJDLFile = open(jdl_file,"w")

tmpCondorJDLFile.write("Universe = vanilla \n")
tmpCondorJDLFile.write("Executable = {} \n".format(executable))
print("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                        .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))
tmpCondorJDLFile.write("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                        .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))

tmpCondorJDLFile.write("Log = {}/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).log \n".format(log_dir,analyzer, sample_name.replace("/","_"), maxjob))
tmpCondorJDLFile.write("Output = {}/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).out \n".format(log_dir, analyzer, sample_name.replace("/","_"), maxjob))
tmpCondorJDLFile.write("Error = {}/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).err \n".format(log_dir, analyzer, sample_name.replace("/","_"), maxjob))

tmpCondorJDLFile.write("+JobQueue=\"Short\" \n")
tmpCondorJDLFile.write("RequestMemory = 4096 \n")
tmpCondorJDLFile.write("RequestCpus = 1 \n")
tmpCondorJDLFile.write("RequestDisk = 4 \n")

tmpCondorJDLFile.write("+RunAsOwner = True \n")
tmpCondorJDLFile.write("+InteractiveUser = true \n")
tmpCondorJDLFile.write("+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9\" \n")
tmpCondorJDLFile.write('+SingularityBindCVMFS = True \n')
tmpCondorJDLFile.write("run_as_owner = True \n")
#tmpCondorJDLFile.write("x509userproxy = {}/x509_proxy \n".format(HOME))
tmpCondorJDLFile.write("should_transfer_files = YES \n")
tmpCondorJDLFile.write("when_to_transfer_output = ON_EXIT_OR_EVICT \n")
tmpCondorJDLFile.write("Transfer_Input_Files = RazorRun,{},bin/Run{},DNNevaluation/EvaluateDNN.py,DNNevaluation/training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5, data/JetVetoMap/{}, data/PileupWeights/{}, data/trigger/{}, data/trigger/{}, Tarballs_for_LPC_Condor/SSLTarball.tar.gz\n".format(inputfilelist,analyzer,jetVetoMap,pileupWeights,HMT,MET))
tmpCondorJDLFile.write("Queue {} \n".format(maxjob))
#tmpCondorJDLFile.write("Queue {} \n".format(2))
tmpCondorJDLFile.close()

if not args.dryRun:
    os.system("condor_submit {} --batch-name {}".format(jdl_file, sample_name))
else:
    print(f"[dryRun] Wrote JDL: {jdl_file}")
    print(f"[dryRun] Would run: condor_submit {jdl_file} --batch-name {sample_name}")
