#version=bparking/V1p19/Data2018_UL/v9/normalized
version=Run3/V1p19/MC_Summer22EE/v1/sixie/v3/normalized/
version=Run3/V1p19/Data2022/v5/normalized/
#ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_1pb_weighted
for sample in \
DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi
do

	input=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/${version}/${sample}.root
	output=/store/user/christiw/displacedJetMuonAnalyzer/${version}
	eval `scram unsetenv -sh`
	gfal-mkdir -p davs://xrootd-redir-stageout.ultralight.org:1095//${output}
	echo "gfal-mkdir -p davs://xrootd-redir-stageout.ultralight.org:1095//${output}"
	gfal-copy -f --checksum-mode=both ${input} davs://xrootd-redir-stageout.ultralight.org:1095/${output}/ 
	echo "gfal-copy -f --checksum-mode=both ${input} davs://xrootd-redir-stageout.ultralight.org:1095/${output}/" 


done
