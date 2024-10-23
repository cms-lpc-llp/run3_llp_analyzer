version=v16
dir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/Data2023/${version}/normalized/
outputRoot=Muon-EXOCSCCluster_Run2023-PromptReco.root


rm -f ${outputRoot}
eval `scram runtime -sh`

hadd ${outputRoot}  ${dir}Muon0-EXOCSCCluster_Run2023B-PromptReco-v1.root  ${dir}Muon1-EXOCSCCluster_Run2023B-PromptReco-v1.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v1.root ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v1.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v2.root  ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v2.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v3.root  ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v3.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v4.root  ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v4.root ${dir}Muon0-EXOCSCCluster_Run2023D-PromptReco-v1.root  ${dir}Muon1-EXOCSCCluster_Run2023D-PromptReco-v1.root ${dir}Muon0-EXOCSCCluster_Run2023D-PromptReco-v2.root  ${dir}Muon1-EXOCSCCluster_Run2023D-PromptReco-v2.root


if [ -f $outputRoot ]
then
        echo "2022 hadd done"
fi

cp ${outputRoot} ${dir}/${outputRoot}

if [ -f $dir/$outputRoot ]
then
        echo "copy succeed"
        rm $outputRoot
fi

dir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/Data2022/${version}/normalized/
outputRoot=DisplacedJet-EXOCSCCluster_Run2022-PromptReco.root
rm -f ${outputRoot}
eval `scram runtime -sh`


hadd ${outputRoot} ${dir}/DisplacedJet-EXOCSCCluster_Run2022E-PromptReco-v1.root ${dir}/DisplacedJet-EXOCSCCluster_Run2022F-PromptReco-v1.root ${dir}/DisplacedJet-EXOCSCCluster_Run2022G-PromptReco-v1.root

if [ -f $outputRoot ]
then
        echo "2023 hadd done"
fi


#eval `scram unsetenv -sh`
#LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
cp ${outputRoot} ${dir}/${outputRoot}

if [ -f $dir/$outputRoot ]
then
	echo "copy succeed"
	rm $outputRoot
fi
