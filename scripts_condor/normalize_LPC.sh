#!/bin/sh

hostname
echo "Job started"
date

isData=$1
sample=$2
inputDir=$3/
outputDir=$4
currentDir=`pwd`
CMSSW_BASE=$5
homeDir=$6
lumi=$7
#user=${homeDir#*/data/}
#user=${homeDir#*/storage/user/}
#runDir=${currentDir}/${user}_${code_dir_suffix}/
runDir=${currentDir}/CMSSW_14_1_0_pre4/src/

normalize_file=llp_${sample}.txt
rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then

	#setup cmssw
	export CWD=${PWD}
	export PATH=${PATH}:/cvmfs/cms.cern.ch/common/
	export SCRAM_ARCH=el8_amd64_gcc12
	echo "PATH: $PATH"
	echo "SCRAM_ARCH: $SCRAM_ARCH"
	scramv1 project CMSSW CMSSW_14_1_0_pre4
	echo "Inside $currentDir:"
	echo "CMSSW_BASE ${CMSSW_BASE}"

	#get grid proxy
	export X509_USER_PROXY=${currentDir}/x509up_u57571
	echo "${currentDir}/x509up_u57571"
	voms-proxy-info
	ls -lah
	cp inputfilelist_${sample}.txt ${runDir} #copy input file list
	cp SSLTarball.tar.gz ${runDir} #copy OpenSSL libraries
	cp xSections.dat ${runDir} #copy xSections.dat file
	cp NormalizeNtuple ${runDir} #copy NormalizeNtuple script
	cp normalization_input.txt ${runDir} #copy normalization input file
	cd ${runDir}
        echo "entering directory: ${runDir}"
	ls
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
	
	#hadd all the jobs for this sample
	#echo "${inputDir}/${sample}*_Job*.root"

	#eosls -f ${inputDir}/ > "inputfilelistForThisJob_${jobnumber}.txt"
	while read line; do xrdcp -f root://cmseos.fnal.gov/${inputDir}/${line} .; done < inputfilelist_${sample}.txt
	#ls ${inputDir}/*_Job*.root | wc
	#echo "${inputDir}/*_Job*.root"
	hadd ${sample}.root *_Job*.root
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
    #eosmkdir -p ${outputDir}
	xrdcp ${runDir}/${output} root://cmseos.fnal.gov/${outputDir}/${output}


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
