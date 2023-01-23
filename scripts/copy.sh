#version=bparking/V1p19/Data2018_UL/v9/normalized
version=Run3/V1p19/Data2022/v1/normalized/
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
