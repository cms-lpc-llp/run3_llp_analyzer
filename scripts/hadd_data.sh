version=v20
dir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/Data2024/${version}/normalized/
outputRoot=Muon-EXOCSCCluster_Run2024-PromptReco.root


rm -f ${outputRoot}
eval `scram runtime -sh`

hadd ${outputRoot} ${dir}Muon0-Run2024B-PromptReco-v1_v1_v1.root  ${dir}Muon0-Run2024E-PromptReco-v2_v1_v1.root ${dir}Muon1-Run2024B-PromptReco-v1_v1_v1.root  ${dir}Muon1-Run2024E-PromptReco-v2_v1_v1.root ${dir}Muon0-Run2024C-PromptReco-v1_v1_v1.root  ${dir}Muon1-Run2024C-PromptReco-v1_v1_v1.root  ${dir}Muon0Run2024F-EXOCSCCluster-PromptReco-v1_v1_v1.root  ${dir}Muon1Run2024F-EXOCSCCluster-PromptReco-v1_v1_v1.root ${dir}Muon0-Run2024D-PromptReco-v1_v1_v1.root  ${dir}Muon1-Run2024D-PromptReco-v1_v1_v1.root  ${dir}Muon0Run2024G-EXOCSCCluster-PromptReco-v1_v1_v1.root  ${dir}Muon1Run2024G-EXOCSCCluster-PromptReco-v1_v1_v1.root ${dir}Muon0-Run2024E-PromptReco-v1_v1_v1.root  ${dir}Muon1-Run2024E-PromptReco-v1_v1_v1.root  


Muon1-Run2024A-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw   21945 Jan 13 18:12 Muon1-Run2024B-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  247781 Jan 13 18:12 Muon1-Run2024C-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  229438 Jan 13 18:12 Muon1-Run2024D-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  175332 Jan 13 18:12 Muon1-Run2024E-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  146453 Jan 13 18:12 Muon1-Run2024E-PromptReco-v2_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  214185 Jan 13 18:12 Muon1-Run2024I-EXOCSCCluster-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  234629 Jan 13 18:12 Muon1-Run2024I-EXOCSCCluster-PromptReco-v2_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  991441 Jan 13 18:12 Muon1Run2024F-EXOCSCCluster-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw 1502380 Jan 13 18:12 Muon1Run2024G-EXOCSCCluster-PromptReco-v1_v1_v1.txt
-rw-rw-r-- 1 christiw christiw  243928 Jan 13 18:12 Muon1Run2024H-EXOCSCCluster-PromptReco-v1_v1_v1.txt

if [ -f $outputRoot ]
then
        echo "2024 hadd done"
fi

cp ${outputRoot} ${dir}/${outputRoot}

if [ -f $dir/$outputRoot ]
then
        echo "copy succeed"
        rm $outputRoot
fi


dir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/Data2023/${version}/normalized/
outputRoot=Muon-EXOCSCCluster_Run2023-PromptReco.root


rm -f ${outputRoot}
eval `scram runtime -sh`

hadd ${outputRoot}  ${dir}Muon0-EXOCSCCluster_Run2023B-PromptReco-v1.root  ${dir}Muon1-EXOCSCCluster_Run2023B-PromptReco-v1.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v1.root ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v1.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v2.root  ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v2.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v3.root  ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v3.root ${dir}Muon0-EXOCSCCluster_Run2023C-PromptReco-v4.root  ${dir}Muon1-EXOCSCCluster_Run2023C-PromptReco-v4.root ${dir}Muon0-EXOCSCCluster_Run2023D-PromptReco-v1.root  ${dir}Muon1-EXOCSCCluster_Run2023D-PromptReco-v1.root ${dir}Muon0-EXOCSCCluster_Run2023D-PromptReco-v2.root  ${dir}Muon1-EXOCSCCluster_Run2023D-PromptReco-v2.root


if [ -f $outputRoot ]
then
        echo "2023 hadd done"
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
