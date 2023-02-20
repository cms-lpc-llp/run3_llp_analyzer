dir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/Data2022/v5/normalized/
outputRoot=DisplacedJet-EXOCSCCluster_Run2022EFG-PromptReco-v1_goodLumi.root

rm -f ${outputRoot}
eval `scram runtime -sh`

hadd ${outputRoot} ${dir}/DisplacedJet-EXOCSCCluster_Run2022E-PromptReco-v1_goodLumi.root ${dir}/DisplacedJet-EXOCSCCluster_Run2022F-PromptReco-v1_goodLumi.root ${dir}/DisplacedJet-EXOCSCCluster_Run2022G-PromptReco-v1_goodLumi.root

if [ -f $outputRoot ]
then
        echo "hadd done"
fi

eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
cp ${outputRoot} ${dir}/${outputRoot}

if [ -f $dir/$outputRoot ]
then
	echo "copy succeed"
	rm $outputRoot
fi
