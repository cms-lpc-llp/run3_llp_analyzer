#!/bin/bash

CACHE=./cache
N_JOBS=128
BIN=../CacheNtuples
while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
        echo "Usage: $0 [options] <input_list_files>"
        echo "Options:"
        echo "  -c, --cache <path>      Path for storing the output caches, defaults to ./cache"
        echo "  -j, --jobs <n>          Number of parallel jobs, defaults to 128"
        echo "  -b, --bin <path>        Path to the binary, defaults to ../CacheNtuples"
        exit 0
        ;;
    -c|--cache)
        CACHE="$2"
        shift
        shift
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
    -|--)
        INP_LIST_FILES="${@:2}"
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

echo "CACHE = $CACHE"
echo "N_JOBS = $N_JOBS"
echo "BIN = $BIN"
echo "INP_LIST_FILES = $INP_LIST_FILES"

if [ ! -f $BIN ]; then
  echo "Binary $BIN not found"
  exit 1
fi

function launch {
    BIN=$1
    FILE=$2
    CACHE=$3

    $BIN $FILE $CACHE > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "$FILE" >> $CACHE/failed.log
        echo "Failed to cache $FILE"
    fi
}

export -f launch

cat $INP_LIST_FILES | parallel -j $N_JOBS --progress $BIN {} $CACHE