#!/bin/bash

BIN=../bin/NormalizeNtuple
TMP_PATH=/tmp/LLP_normalize_tmpfiles
N_JOBS="$nproc"
LUMI=1
RSYNC=

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
        echo "Usage: $0 [options] -- [paths containing files to be combined & normalized]"
        echo "Options:"
        echo "  -j, --jobs <number>     Number of parallel jobs, default $nproc"
        echo "  -b, --bin <path>        Path to the binary, default ../bin/NormalizeNtuple"
        echo "  -l, --lumi <number>     Integrated luminosity, default 1 (pb)"
        exit 0
        ;;
    -j|--jobs)
        N_JOBS="$2"
        shift
        shift
        ;;
    -b|--bin)
        BIN="$2"
        shift
        shift
        ;;
    -l|--lumi)
        LUMI="$2"
        shift
        shift
        ;;
    -t|--temp)
        TMP_PATH="$2"
        shift
        shift
        ;;
    -|--)
        INP_PATHS="${@:2}"
        break
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

echo "N_JOBS: $N_JOBS"
echo "BIN: $BIN"
echo "LUMI: $LUMI pb"
echo "TMP_PATH: $TMP_PATH"
echo "INP_PATHS: $INP_PATHS"

for INP_PATH in $INP_PATHS; do
    if [ ! -d "$INP_PATH" ]; then
        echo INP_PATH $INP_PATH is not a directory, exiting
        exit 1
    fi
done

if [ ! -f $BIN ]; then
    echo "Binary $BIN not found, exiting"
    exit 1
fi

if [ ! -f "$(dirname $BIN)/data/xSections.dat" ]; then
    echo "xSections.dat not found at $(dirname $BIN)/data/xSections.dat, exiting"
    exit 1
fi

function launch {

    INP_PATH=$(realpath $1)
    BIN=$2
    LUMI=$3
    IS_DATA=$(echo $INP_PATH | grep -q "Data" && echo "yes" || echo "no")
    cd $INP_PATH
    if [ -f $INP_PATH.root ]; then
        rm $INP_PATH.root
    fi
    hadd $INP_PATH.root $(ls *.root | sort -h | tr '\n' ' ')

    if [ $? -ne 0 ]; then
        echo "Hadd failed"
        exit 1
    fi

    if [ $IS_DATA == "no" ]; then
        cd $(dirname $BIN)
        TMP_FILE=$TMP_PATH/$(basename $INP_PATH).txt
        echo "$(basename $INP_PATH) $INP_PATH.root" > $TMP_FILE
        echo $BIN $TMP_FILE $LUMI
        if [ $? -ne 0 ]; then
            echo "Normalization failed"
            exit 1
        fi
        rm $TMP_FILE
        rm $INP_PATH.root
        # rm -rf $INP_PATH
    fi
}

export -f launch

echo $INP_PATHS | parallel --progress -j $N_JOBS launch {} $BIN $LUMI
