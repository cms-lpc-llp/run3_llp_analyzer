#!/bin/bash

N_JOBS="$(nproc)"
_CMSSW_BASE=$(realpath ~/repo/LLP/CMSSW_9_4_4)
JSON_DIR=../data/Certification/
BIN=

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
        echo "Usage: $0 [options]"
        echo "Options:"
        echo "  -j, --jobs <number>     Number of parallel jobs, default $\(nproc\)"
        echo "  -e, --env  <cmssw path to be used> Path to the CMSSW base, default ~/repo/LLP/CMSSW_9_4_4"
        echo "  -d, --data <path>       Path to the data files"
        echo "  -b, --bin  <path>       Path to the binary to be used"
        exit 0
        ;;
    -j|--jobs)
        N_JOBS="$2"
        shift
        shift
        ;;
    -d|--data)
        DATA_PATH="$2"
        shift
        shift
        ;;
    -e|--env)
        _CMSSW_BASE="$2"
        shift
        shift
        ;;
    -jd|--json-dir)
        JSON_DIR="$2"
        shift
        shift
        ;;
    -b|--bin)
        BIN=$(realpath "$2")
        shift
        shift
        ;;
    -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;
    *)
        echo "Unknown option $1"
        exit 1
        ;;
  esac
done


function write_loader {

    JSON_PATH=$(realpath $1)
    OUT_PATH=$(mktemp --suffix=.py)

    PYTHON="import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

# setup process
process = cms.Process('FWLitePlots')
process.inputs = cms.PSet (
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
)

# get JSON file correctly parced
JSONfile = '${JSON_PATH}'
myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

process.inputs.lumisToProcess.extend(myList)"

    echo "$PYTHON" > $OUT_PATH

    echo $OUT_PATH
}

function get_json_path {
    # This function should be adapted to the actual path of the golden JSON files
    INP_FILE=$(realpath $1)
    JSON_DIR=$2
    YEAR=$(echo $INP_FILE | grep -o "20[0-9]\{2\}[A-Z]\?")
    YEAR=${YEAR:0:4}

    ERA=${RUN:7:1}

    JSON_NAME=Cert_Collisions${YEAR}_*_Golden.json

    JSON_PATH=$(find $JSON_DIR -name $JSON_NAME)

    if [ -z $JSON_PATH ]; then
        echo "******************************************************************" 1>&2
        echo "JSON path not found for file $INP_FILE:" 1>&2
        echo "  CMD: find $JSON_DIR -name $JSON_NAME" 1>&2
        exit 1
    fi

    echo $JSON_PATH
}

function filter {
    INP_FILE=$(realpath $1)
    JSON_FILE=$(get_json_path $INP_FILE $2)
    BIN=$3
    OUT_FILE=${INP_FILE%.root*}_goodlumi.root

    LOADER=$(write_loader $JSON_FILE)
    $BIN $LOADER $INP_FILE $OUT_FILE
    rm $LOADER
}

export -f write_loader
export -f get_json_path
export -f filter

if [ -z "$CMSSW_VERSION" ] || [ $CMSSW_VERSION != "CMSSW_9_4_4" ]; then
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd "$_CMSSW_BASE"
    cmsenv
    if [ $? -ne 0 ]; then
        echo "CMSSW_9_9_4 is not activated, and it does not exist at specified path $_CMSSW_BASE, exiting"
        exit 1
    fi
    cd -
fi

if [ -z "$CMSSW_VERSION" ] || [ $CMSSW_VERSION != "CMSSW_9_4_4" ]; then
    echo "CMSSW_BASE not set to CMSSW_9_4_4, exiting"
    exit 1
fi

if [ ! -d $DATA_PATH ]; then
    echo "Data path $DATA_PATH not found, exiting"
    exit 1
fi

ls $DATA_PATH/*.root | grep "20[0-9]\{2\}[A-Z]\?" | parallel --progress -j $N_JOBS filter {} $JSON_DIR $BIN
