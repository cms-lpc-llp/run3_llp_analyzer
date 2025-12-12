#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`
analysisType=$1
inputfilelist=$2
isData=$3
filePerJob=$4
jobnumber=$5
maxjob=$6
sample=${inputfilelist##*/}
sample=${sample%.txt}
outputfile=${sample}_Job${jobnumber}_of_${maxjob}.root
outputDirectory=$7
analyzerTag=$8
CMSSW_BASE=$9
homeDir=${10}
option=${11}
currentDir=`pwd`
echo "currentDir: ${currentDir}"
#user=${homeDir#*/data/}
#user=${homeDir#*/storage/user/}
#runDir=${currentDir}/${analyzerTag}/
runDir=${currentDir}/CMSSW_14_1_0_pre4/src/

rm -rf ${runDir}
#mkdir -p ${runDir}
echo ${CMSSW_BASE}
echo homeDir: ${homeDir}

:'
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	#setup cmssw
	ls -la
	cd CMSSW_10_6_20/src/
	#workDir=`pwd`
	#echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	eval `scram runtime -sh`
	cd ${currentDir}
	#echo `which root`

	#cd ${runDir}
	#echo "entering directory: ${runDir}"
	#echo "${CMSSW_BASE}/src/run3_llp_analyzer/RazorRun"
'
#copying LPC commands

echo ${CMSSW_BASE}
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	export CWD=${PWD}
	export PATH=${PATH}:/cvmfs/cms.cern.ch/common/
	export SCRAM_ARCH=el8_amd64_gcc12
	echo "PATH: $PATH"
	echo "SCRAM_ARCH: $SCRAM_ARCH"
	scramv1 project CMSSW CMSSW_14_1_0_pre4
	#cmsrel CMSSW_10_6_20
	echo "Inside $currentDir:"
	ls -lah
	#cp RazorRun CMSSW_10_6_20/src/run3_llp_analyzer/
	echo "CMSSW_BASE ${CMSSW_BASE}"
	#cd ${CMSSW_BASE}
	#cd CMSSW_10_6_20/src
	ls -lah

#cp RazorRun CMSSW_10_6_20/src
#cp Runllp_MuonSystem_CA_TnP CMSSW_10_6_20/src
#cp EvaluateDNN.py CMSSW_10_6_20/src
#cp training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5
	if [ -f *"Runllp"* ]
	then
		#cp $CMSSW_BASE/src/run3_llp_analyzer/RazorRun ./
		#cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/EvaluateDNN.py ./
        #        cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/*.h5 .
		cp Run${analysisType} ${runDir}
		cp EvaluateDNN.py ${runDir}
		cp training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5 ${runDir}
		cp ${sample}.txt ${runDir}
		cp *Summer*.root ${runDir} #should get Pileupreweight files
		cp *Prompt*.root ${runDir} #should get jet veto files
		cp *fficiencies*.root ${runDir} #should get HMT trigger efficiencies files
		cp *MET*.root ${runDir} #should get MET trigger efficiency files
		cp SSLTarball.tar.gz ${runDir} #copy OpenSSL libraries

		


		#get grid proxy
		export X509_USER_PROXY=${currentDir}/x509up_u57571
		echo "${currentDir}/x509up_u57571"
		voms-proxy-info

		cd ${runDir}
		tar -xvf SSLTarball.tar.gz
		export OPENSSL_HOME=${runDir}/openssl-1.1
		export LD_LIBRARY_PATH=${OPENSSL_HOME}/lib:${LD_LIBRARY_PATH}
		export PATH=${OPENSSL_HOME}/bin:${PATH}
		echo "PATH: $PATH"
		echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
		echo "entering directory: ${runDir}"
		#cmsenv
		eval `scram runtime -sh`
		source /cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.30.07-024df6516c17fd2edef848a927a788f1/bin/thisroot.sh
		#run the job
		# echo "cat ${inputfilelist} | awk \"NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})\" > inputfilelistForThisJob_${jobnumber}.txt"
		cat ${sample}.txt | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
                echo "Running on these input files:"
                cat inputfilelistForThisJob_${jobnumber}.txt
                echo "************************************"
	
		echo "start copying files"
		while read line; do xrdcp -f ${line/xrootd/eos} .; done < "inputfilelistForThisJob_${jobnumber}.txt"
		rm inputfilelistForThisJob_${jobnumber}.txt
		#ls -ltr *ntupler*.root
		ls -ltr
		#renaming files with -- in the name
		#for filename in *; do
		#	if [ ${filename:0:2}=="--" ]; then
		#		mv $filename ${filename:2}
		#	fi
		#done
		echo "now sending ls output to new file"
		#ls *NANO*.root
		ls ./*NANO*.root ./*Job*.root ./*EXO*.root > inputfilelistForThisJob_${jobnumber}.txt
		
		
		

		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""
		echo " "; echo "Starting razor run job now"; echo " ";
		if [ ${analysisType} == "MakeMCPileupDistribution" ]
		then
			echo "./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -f=${outputfile}"
			./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -f=${outputfile}
		else
			#echo ./Runllp_MuonSystem_CA_TnP inputfilelistForThisJob_${jobnumber}.txt -d=${isData}  -f=${outputfile} -l=${analyzerTag}

			if [[ ${isData} == "no" ]]; then #check if it is DY MC, if a signal sample this won't work
				echo "Running on MC sample"
				echo ./Run${analysisType} inputfilelistForThisJob_${jobnumber}.txt  -f=${outputfile} -l=${analyzerTag} -n=${option}
				./Run${analysisType} inputfilelistForThisJob_${jobnumber}.txt  -f=${outputfile} -l=${analyzerTag} -n=${option}
			else
				echo "Running on data sample"
				echo ./Run${analysisType} inputfilelistForThisJob_${jobnumber}.txt  --isData -f=${outputfile} -l=${analyzerTag} -n=${option}
				./Run${analysisType} inputfilelistForThisJob_${jobnumber}.txt  --isData -f=${outputfile} -l=${analyzerTag} -n=${option}
			fi
			# #./Runllp_MuonSystem_CA_TnP inputfilelistForThisJob_${jobnumber}.txt  --isData  -f=${outputfile}
			# echo ./${analysisType} inputfilelistForThisJob_${jobnumber}.txt  -f=${outputfile} -l=${analyzerTag}
			# ./${analysisType} inputfilelistForThisJob_${jobnumber}.txt  -f=${outputfile} -l=${analyzerTag}
		fi
			
		echo ${outputfile}
		echo ${outputDirectory}
		ls *Job*.root > output.txt
		echo "Output ROOT files: "
		cat output.txt
		##^_^##
		echo "RazorRun_T2 finished"
		date

		#COMMENTED OUT FOR PROMPT HNL - NO CLUSTERS
		if [[ ${analysisType} != *"noClusters"* && ${analysisType} != *"TrigEff"* ]]; then
			echo "start DNN evaluation"
			source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
			python EvaluateDNN.py --in_file ${outputfile}
		fi
		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		#mkdir -p ${outputDirectory}
		echo "trying to copy file"
		xrdcp -f ${outputfile} root://cmseos.fnal.gov/${outputDirectory}/${outputfile}
		echo "copied file"
		:'
		if [ -f ${outputDirectory}/${outputfile} ]
		then
			echo ${outputfile} "copied"
		else
			echo ${outputfile} "not copied"
		fi
		'
	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
#rm -rf ${runDir}
echo "Job finished"
date
end_time=`date +%s`
runtime=$((end_time-start_time))
echo ${runtime}
