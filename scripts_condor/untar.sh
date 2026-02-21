#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 11 ]; then
  echo "Usage: $0 <analysisType> <inputfilelist> <isData> <filePerJob> <jobnumber> <maxjob> <outputDirectory> <analyzerTag> <cmsswBaseHint> <homeDir> <cmsswTarball>" >&2
  exit 2
fi

analysisType=$1
inputfilelist=$2
isData=$3
filePerJob=$4
jobnumber=$5
maxjob=$6
outputDirectory=$7
analyzerTag=$8
cmsswBaseHint=$9
homeDir=${10}
cmsswTarball=${11}

echo "wrapper pwd=$(pwd)"
echo "analysisType=${analysisType}"
echo "inputfilelist=${inputfilelist}"
echo "submitHostCmsswHint=${cmsswBaseHint} (informational only)"
echo "cmsswTarball=${cmsswTarball}"

# Condor may pass a path; use basename fallback if needed.
if [ ! -f "${cmsswTarball}" ]; then
  cmsswTarball="$(basename "${cmsswTarball}")"
fi
if [ ! -f "${cmsswTarball}" ]; then
  echo "CMSSW tarball not found: ${cmsswTarball}" >&2
  ls -la
  exit 2
fi

if [ ! -f "runAnalyzer_ucsd.sh" ]; then
  echo "runAnalyzer_ucsd.sh not found in job sandbox" >&2
  ls -la
  exit 2
fi

echo "Extracting CMSSW tarball: ${cmsswTarball}"
top_dir="$(tar -tzf "${cmsswTarball}" | sed -n '1s#/.*##p')"
if [ -z "${top_dir}" ]; then
  echo "Failed to detect top-level directory in tarball ${cmsswTarball}" >&2
  exit 2
fi

tar -xzf "${cmsswTarball}"

local_cmssw_base="${PWD}/${top_dir}"
if [ ! -d "${local_cmssw_base}/src" ]; then
  echo "Extracted CMSSW layout invalid: ${local_cmssw_base}" >&2
  ls -la
  exit 2
fi

echo "Using extracted CMSSW_BASE=${local_cmssw_base}"

exec bash runAnalyzer_ucsd.sh \
  "${analysisType}" \
  "${inputfilelist}" \
  "${isData}" \
  "${filePerJob}" \
  "${jobnumber}" \
  "${maxjob}" \
  "${outputDirectory}" \
  "${analyzerTag}" \
  "${local_cmssw_base}" \
  "${homeDir}"
