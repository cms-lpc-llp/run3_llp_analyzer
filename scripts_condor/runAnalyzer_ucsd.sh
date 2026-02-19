#!/bin/bash

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
RUN_DIAG=${RUN_DIAG:-0}

homeDirNoSlash=${homeDir%/}
user=`basename "${homeDirNoSlash}"`
runDir=${currentDir}/${user}_${analyzerTag}/

rm -rf ${runDir}
mkdir -p ${runDir}
echo ${CMSSW_BASE}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	# setup cmssw
	cd ${CMSSW_BASE}/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	# TODO(uaf): confirm SCRAM_ARCH for UAF (el8/el9) if needed.
	export SCRAM_ARCH=el8_amd64_gcc12
	if ! eval `scramv1 runtime -sh`; then
		echo "scram runtime failed" >&2
		exit 2
	fi

	# Debug: environment snapshot (stdout + stderr)
	echo "=== system ==="
	uname -a || true
	cat /etc/os-release || true
	echo "=== env ==="
	echo "CMSSW_BASE=${CMSSW_BASE}"
	echo "CMSSW_RELEASE_BASE=${CMSSW_RELEASE_BASE}"
	echo "SCRAM_ARCH=${SCRAM_ARCH}"
	echo "PATH=${PATH}"
	echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

	# Ensure CMSSW bins/libs are visible
	if [ -z "${CMSSW_RELEASE_BASE:-}" ]; then
		export CMSSW_RELEASE_BASE="/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_14_1_0_pre4"
		echo "CMSSW_RELEASE_BASE set to ${CMSSW_RELEASE_BASE}"
	fi

	export PATH="${CMSSW_BASE}/bin/${SCRAM_ARCH}:${CMSSW_BASE}/external/${SCRAM_ARCH}/bin:${CMSSW_RELEASE_BASE}/bin/${SCRAM_ARCH}:${CMSSW_RELEASE_BASE}/external/${SCRAM_ARCH}/bin:${PATH}"
	export LD_LIBRARY_PATH="${CMSSW_BASE}/external/${SCRAM_ARCH}/lib:${CMSSW_BASE}/external/${SCRAM_ARCH}/lib64:${CMSSW_BASE}/lib/${SCRAM_ARCH}:${CMSSW_RELEASE_BASE}/external/${SCRAM_ARCH}/lib:${CMSSW_RELEASE_BASE}/external/${SCRAM_ARCH}/lib64:${CMSSW_RELEASE_BASE}/lib/${SCRAM_ARCH}:${LD_LIBRARY_PATH}"

	echo "root=$(command -v root)"

	cd ${runDir}
	echo "entering directory: ${runDir}"
	echo "${CMSSW_BASE}/src/run3_llp_analyzer/RazorRun"
	if [ -f ${CMSSW_BASE}/src/run3_llp_analyzer/RazorRun ]
	then
		cp $CMSSW_BASE/src/run3_llp_analyzer/RazorRun ./
		# cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/EvaluateDNN.py ./
		# cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/*.h5 .

		# get grid proxy
		if [ -n "${X509_USER_PROXY}" ]; then
			echo "Using X509_USER_PROXY=${X509_USER_PROXY}"
		elif [ -f "/tmp/x509up_u$(id -u)" ]; then
			export X509_USER_PROXY="/tmp/x509up_u$(id -u)"
			echo "Using X509_USER_PROXY=${X509_USER_PROXY}"
		else
			export X509_USER_PROXY="${homeDirNoSlash}/x509_proxy"
			echo "Using X509_USER_PROXY=${X509_USER_PROXY}"
		fi
		voms-proxy-info

		# Stage required helper ROOT files into the run directory
		echo "Staging helper ROOT files into ${runDir}"
		helper_files=(
			"$CMSSW_BASE/src/run3_llp_analyzer/data/trigger/METTriggerEff_Summer24.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/trigger/HMT_Efficiencies_2024.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/pileup/PileupReweight_Summer24.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/pileup/PileupReweight_Summer23BPix.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/pileup/PileupReweight_Summer23.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/pileup/PileupReweight_Summer22EE.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/pileup/PileupReweight_Summer22.root"
			"$CMSSW_BASE/src/run3_llp_analyzer/data/jetvetomaps/Summer24Prompt24_2024BCDEFGHI.root"
		)
		for f in "${helper_files[@]}"; do
			if [ -f "$f" ]; then
				cp -f "$f" "${runDir}/"
			else
				echo "missing helper file: $f" >&2
			fi
		done

		# run the job
		cat ${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo "=== input list diagnostics ==="
		echo "inputfilelist=${inputfilelist}"
		echo "jobnumber=${jobnumber} filePerJob=${filePerJob}"
		echo "first 3 files in list:"
		head -n 3 inputfilelistForThisJob_${jobnumber}.txt || true
		first_file=$(head -n 1 inputfilelistForThisJob_${jobnumber}.txt || true)
		echo "first_file=${first_file}"

		if [ "${RUN_DIAG}" = "1" ]; then
			echo "=== runtime diagnostics (RUN_DIAG=1) ==="
			echo "pwd=$(pwd)"
			echo "whoami=$(whoami)"
			echo "id=$(id)"
			echo "ls -la runDir:"
			ls -la
			echo "ls -la CMSSW bin:"
			ls -la "${CMSSW_BASE}/bin/${SCRAM_ARCH}" || true
			echo "ls -la CMSSW lib:"
			ls -la "${CMSSW_BASE}/lib/${SCRAM_ARCH}" || true
		fi

		# ROOT inspection of the first input file (helps catch schema mismatch)
		if [ -n "${first_file}" ] && command -v root >/dev/null 2>&1; then
			echo "=== ROOT file inspection (first file) ==="
			root -l -b -q -e "
auto f = std::unique_ptr<TFile>{TFile::Open(\"${first_file}\")};
if (!f || f->IsZombie()) { std::cerr << \"Failed to open file\" << std::endl; gSystem->Exit(2); }
std::cout << \"Keys:\" << std::endl;
f->GetListOfKeys()->Print();
TTree* t = (TTree*)f->Get(\"MuonSystem\");
if (t) {
  std::cout << \"MuonSystem entries=\" << t->GetEntries() << std::endl;
  TBranch* b1 = t->GetBranch(\"nTaus\");
  TBranch* b2 = t->GetBranch(\"tauE\");
  if (b1) { b1->GetListOfLeaves()->Print(); } else { std::cout << \"Missing branch nTaus\" << std::endl; }
  if (b2) { b2->GetListOfLeaves()->Print(); } else { std::cout << \"Missing branch tauE\" << std::endl; }
  std::cout << \"First 25 branches:\" << std::endl;
  auto br = t->GetListOfBranches();
  int n = br ? br->GetEntries() : 0;
  for (int i=0; i<std::min(25, n); ++i) {
    std::cout << \"  \" << br->At(i)->GetName() << std::endl;
  }
} else {
  std::cout << \"MuonSystem tree not found\" << std::endl;
}
";
		else
			echo "root not available or first file missing; skipping ROOT inspection"
		fi
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"

		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""
		echo " "; echo "Starting razor run job now"; echo " ";

		cp $CMSSW_BASE/src/run3_llp_analyzer/bin/Runllp_MuonSystem_CA_mdsnano ./
		if [ ! -x "./Runllp_MuonSystem_CA_mdsnano" ]; then
			echo "Runllp binary missing or not executable" >&2
			ls -l $CMSSW_BASE/src/run3_llp_analyzer/bin/Runllp_MuonSystem_CA_mdsnano || true
			exit 3
		fi
		if command -v ldd >/dev/null 2>&1; then
			ldd ./Runllp_MuonSystem_CA_mdsnano | egrep "not found|tbb" >&2 || true
		fi

		if [ "${analysisType}" = "MakeMCPileupDistribution" ]
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

		# echo "start DNN evaluation"
		# # TODO(uaf): confirm LCG view for UAF if needed.
		# source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh
		# python EvaluateDNN.py --in_file ${outputfile}

		sleep 2
		echo "I slept for 2 second"

		## job finished, copy file to output location
		echo "copying output file to ${outputDirectory}"
		if [ ! -f "${outputfile}" ]; then
			echo "output file missing: ${outputfile}" >&2
			ls -l
			exit 4
		fi

		# Resolve output destination
		if [[ "${outputDirectory}" == root://* ]]; then
			dest="${outputDirectory%/}/${outputfile}"
		elif [[ "${outputDirectory}" == /ceph/cms/* ]]; then
			dest="root://redirector.t2.ucsd.edu:1095/${outputDirectory#/ceph/cms}/${outputfile}"
		elif [[ "${outputDirectory}" == /store/* ]]; then
			dest="root://redirector.t2.ucsd.edu:1095${outputDirectory%/}/${outputfile}"
		else
			dest="${outputDirectory%/}/${outputfile}"
		fi

		echo "copying output file via gfal-copy to ${dest}"
		if command -v gfal-copy >/dev/null 2>&1; then
			gfal-copy -f -p -t 1800 "file://${runDir}/${outputfile}" "${dest}"
			gfal_rc=$?
		else
			gfal_rc=127
		fi

		if [ ${gfal_rc} -ne 0 ]; then
			echo "gfal-copy failed (rc=${gfal_rc}), trying xrdcp" >&2
			if command -v xrdcp >/dev/null 2>&1; then
				xrdcp -f "${outputfile}" "${dest}"
				xrd_rc=$?
			else
				xrd_rc=127
			fi
			if [ ${xrd_rc} -ne 0 ]; then
				echo ${outputfile} "not copied" >&2
			else
				echo ${outputfile} "copied"
			fi
		else
			echo ${outputfile} "copied"
		fi

		# Clean CMSSW env after copy
		eval `scram unsetenv -sh`



	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
echo "Job finished"
date
end_time=`date +%s`
runtime=$((end_time-start_time))
echo ${runtime}
