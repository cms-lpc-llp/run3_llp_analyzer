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
currentDir=`pwd`
#user=${homeDir#*/data/}
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_${analyzerTag}/

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
	echo "${CMSSW_BASE}/src/run3_llp_analyzer/RazorRun"
	if [ -f ${CMSSW_BASE}/src/run3_llp_analyzer/RazorRun ]
	then
		cp $CMSSW_BASE/src/run3_llp_analyzer/RazorRun ./
		cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/EvaluateDNN.py ./
                cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/*.h5 .
		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy
		echo "${homeDir}x509_proxy"
		voms-proxy-info


		#run the job
		# echo "cat ${inputfilelist} | awk \"NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})\" > inputfilelistForThisJob_${jobnumber}.txt"
		cat ${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
                echo "Running on these input files:"
                cat inputfilelistForThisJob_${jobnumber}.txt
                echo "************************************"
	
		echo "start copying files"
		while read line; do   xrdcp ${line} .; done < "inputfilelistForThisJob_${jobnumber}.txt"
		rm inputfilelistForThisJob_${jobnumber}.txt
		ls displaced*.root > inputfilelistForThisJob_${jobnumber}.txt
		
		

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
			echo ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData}  -f=${outputfile} -l=${analyzerTag}
			./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData}  -f=${outputfile} -l=${analyzerTag}
		fi

		echo ${outputfile}
		echo ${outputDirectory}
		ls *root > output.txt
		echo "Output ROOT files: "
		cat output.txt
		##^_^##
		echo "RazorRun_T2 finished"
		date

		echo "start DNN evaluation"
		source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh
		python EvaluateDNN.py --in_file ${outputfile}
		
		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		mkdir -p ${outputDirectory}
		cp ${outputfile} ${outputDirectory}/${outputfile}
		if [ -f ${outputDirectory}/${outputfile} ]
		then
			echo ${outputfile} "copied"
		else
			echo ${outputfile} "not copied"
		fi

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
