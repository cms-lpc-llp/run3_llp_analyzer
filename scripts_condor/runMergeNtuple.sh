#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`
inputNtupleList=$1
ntuplePerJob=$2
ntuple_index=$3
maxNtupleJob=$4
inputNanoList=$5
nanoPerJob=$6
nanoaod_index=$7
maxNanoJob=$8
total_jobnumber=$9
sample=${inputNtupleList##*/}
sample=${sample%.txt}

maxjob=`python -c "print int(${maxNanoJob}*${maxNtupleJob})"`
outputfile=${sample}_Job${total_jobnumber}_of_${maxjob}_NtupleJob${ntuple_index}_Of_${maxNtupleJob}_NanoJob${nanoaod_index}_Of_${maxNanoJob}.root
outputDirectory=${10}
CMSSW_BASE=${11}
homeDir=${12}
currentDir=`pwd`
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_MergeNtuple/


rm -rf ${runDir}
mkdir -p ${runDir}
echo ${CMSSW_BASE}
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	#setup cmssw
	cd ${CMSSW_BASE}/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	echo "${CMSSW_BASE}/src/run3_llp_analyzer/MergeNtuples"
	if [ -f ${CMSSW_BASE}/src/run3_llp_analyzer/MergeNtuples ]
	then
		cp $CMSSW_BASE/src/run3_llp_analyzer/MergeNtuples ./
		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy
		echo "${homeDir}x509_proxy"
		voms-proxy-info



		#run the job
		cat ${inputNtupleList} | awk "NR > (${ntuple_index}*${ntuplePerJob}) && NR <= ((${ntuple_index}+1)*${ntuplePerJob})" > inputNtupleList.txt
		cat ${inputNanoList} | awk "NR > (${nanoaod_index}*${nanoPerJob}) && NR <= ((${nanoaod_index}+1)*${nanoPerJob})" > inputNanoList.txt
		echo ""
		echo "************************************"
		echo "Running on these input NTUPLE files:"
		cat inputNtupleList.txt
		echo "************************************"
		echo ""

		echo ""
		echo "************************************"
		echo "Running on these input NANOAOD files:"
		cat inputNanoList.txt
		echo "************************************"
		echo ""


		echo " "; echo "Starting MergeNtuples job now"; echo " ";
		echo ./MergeNtuples inputNtupleList.txt inputNanoList.txt ${outputfile} ""
		./MergeNtuples inputNtupleList.txt inputNanoList.txt ${outputfile} ""
		

		echo ${outputfile}
		echo ${outputDirectory}
		ls *root > output.txt
		echo "Output ROOT files: "
		cat output.txt
		##^_^##
		echo "MergeNtuple finished"
		date

		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		mkdir -p ${outputDirectory}
		while IFS= read -r line
		do
        		echo $line
			cp ${line} ${outputDirectory}/${outputfile}
			if [ -f ${outputDirectory}/${line} ]
			then
				echo ${line} "copied"
			else
				echo ${line} "not copied"
			fi
		done <"output.txt"

	else
		echo echo "WWWWYYYY ============= failed to access file MergeNtuple, job anandoned"
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
