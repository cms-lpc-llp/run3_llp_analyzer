#!/bin/sh

hostname
echo "Job started"
date

isData=$1
sample=$2
inputDir=$3/${sample}/
outputDir=$4
currentDir=`pwd`
CMSSW_BASE=$5
homeDir=$6
lumi=$7
#user=${homeDir#*/data/}
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_${code_dir_suffix}/

normalize_file=llp_${sample}.txt
rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then

	#setup cmssw
	cd ${CMSSW_BASE}/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	ulimit -c 0
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
        echo "entering directory: ${runDir}"

	#hadd all the jobs for this sample
	echo "${inputDir}/${sample}*_Job*.root"
	ls ${inputDir}/*_Job*.root | wc
	echo "${inputDir}/*_Job*.root"
	hadd ${sample}.root ${inputDir}/*_Job*.root
	output=${sample}.root
        ls

	echo "start normalization"
	if [ ${isData} == "no" ]
        then
		eval `scramv1 runtime -sh`
		echo "$CMSSW_BASE/src/run3_llp_analyzer/data/xSections.dat"
		if [ -f $CMSSW_BASE/src/run3_llp_analyzer/data/xSections.dat ]
		then
			mkdir -p data
			cp $CMSSW_BASE/src/run3_llp_analyzer/data/xSections.dat data/xSections.dat
		else
			echo "data/xSections.dat doesn't exist"

		fi

		#create normalization file
		rm -f $normalize_file
		echo "${sample} ${runDir}/${output}" > $normalize_file
		cat $normalize_file

		if [ -f $normalize_file ]
		then
			echo "normalization file created"
		fi

		#normalize
		if [ -f $CMSSW_BASE/src/run3_llp_analyzer/NormalizeNtuple ]
        	then
        	        cp $CMSSW_BASE/src/run3_llp_analyzer/NormalizeNtuple ./
		        ./NormalizeNtuple ${normalize_file} ${lumi}
		else
			echo "NormalizeNtuple not found"
		fi
		echo "Normalization done"
		ls
		output=${output%.root*}_${lumi}pb_weighted.root
		echo "Output file name: ${output}"
	fi
	sleep 2
        echo "I slept for 2 second"



	# copy normalized file back to hadoop
        mkdir -p ${outputDir}
	cp ${runDir}/${output} ${outputDir}/${output}


	if [ -f ${outputDir}/${output} ]
	then
	        echo "copied succeed"
	else
	        echo "copied failed"
	fi


else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
#rm -rf ${runDir}
echo "Job finished"
date
