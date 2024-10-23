#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`


inputfilelist=$1
filePerJob=$2
jobnumber=$3
maxjob=$4
sample=${inputfilelist##*/}
sample=${sample%.txt}
outputfile=${sample}_Job${jobnumber}_of_${maxjob}
outputDirectory=$5
CMSSW_BASE=$6
homeDir=$7
currentDir=`pwd`
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
	echo "${CMSSW_BASE}/src/run3_llp_analyzer/MergeNtuples"
	if [ -f ${CMSSW_BASE}/src/run3_llp_analyzer/MergeNtuples ]
	then
		cp $CMSSW_BASE/src/run3_llp_analyzer/MergeNtuples ./
		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy
		echo "${homeDir}x509_proxy"
		voms-proxy-info
		


		#run the job
		cat ${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		
		echo ""
		echo "************************************"
		echo "Running on these commands:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""

		list_base_path=${CMSSW_BASE}/src/run3_llp_analyzer/scripts_local/

		test=`cat inputfilelistForThisJob_${jobnumber}.txt`
		stringarray=($test)
		inputNtupleList=${list_base_path}/${stringarray[0]}
		inputNanoList=${list_base_path}/${stringarray[1]}
		
		
		echo ${inputNtupleList}
		echo "************************************"
		echo "Running on these input NTUPLE files:"
		cat ${inputNtupleList}
		echo "************************************"
		echo ""

		echo "start copying files"

		while read line || [ -n "$line" ]; do echo $line; done < ${inputNtupleList}
		while read line || [ -n "$line" ]; do xrdcp ${line} .; done < ${inputNtupleList}
        ls *.root > inputNtupleList.txt
		inputNtupleList=inputNtupleList.txt
		
		echo ""
		echo "************************************"
		echo "Running on these input NTUPLE files:"
		cat ${inputNtupleList}
		echo "************************************"
		echo ""

		echo ""
		echo "************************************"
		echo "Running on these input NANOAOD files:"
		cat ${inputNanoList}
		echo "************************************"
		echo ""
		split -d -l 1 ${inputNtupleList}  tempNtuple_ --additional-suffix=.txt --suffix-length=3 --numeric-suffixes=1

		# Nlines=$(wc -l < "${inputNtupleList}")
		Nlines=$(ls "tempNtuple_"*txt | wc -l)

		echo "temp ntuple files, ${Nlines} in ntuple list";
		ls

		echo " ";echo "Starting MergeNtuples job loop"; echo " ";

		for i in $( seq -f "%03g" 1 ${Nlines} )
		do
			echo " "; echo "Starting MergeNtuples job now for Line ${i}"; echo " ";
			cat tempNtuple_${i}.txt
			echo "************************************"

			echo ./MergeNtuples tempNtuple_${i}.txt ${inputNanoList} ${outputfile}_${i}.root "";
			./MergeNtuples tempNtuple_${i}.txt ${inputNanoList} ${outputfile}_${i}.root "";
		done

		

		echo ${outputfile}
		echo ${outputDirectory}
		ls ${outputfile}_*root > output.txt
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
			cp ${line} ${outputDirectory}/${line}
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
