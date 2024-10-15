#!/bin/bash

LINE=$1
ARGS_FILE=$2
BIN=$3

ARGS=$(sed -n "${LINE}p" $ARGS_FILE)
INP_NTUPLE_LIST=$(echo $ARGS | awk '{print $1}')
INP_NANO_LIST=$(echo $ARGS | awk '{print $2}')
OUT_PATH=$(echo $ARGS | awk '{print $3}')
CACHE_PATH=$(echo $ARGS | awk '{print $4}')

LOG_PATH=${OUT_PATH%.root}.log

if [ -f $OUT_PATH ] || [ -f $OUT_PATH.nomatch ]; then
    exit 0
fi

mkdir -p $(dirname $OUT_PATH)

echo 'start merging' > $LOG_PATH
echo "TMP_PATH=$TMP" >> $LOG_PATH

$BIN $INP_NTUPLE_LIST "$INP_NANO_LIST" $OUT_PATH $CACHE_PATH >> $LOG_PATH 2>&1

if [ $? -ne 0 ] || [ ! -f $OUT_PATH ]; then
    # rm -r $TMP
    echo "Failed to merge $INP_NTUPLE_LIST and $INP_NANO_LIST" 1>&2
    echo "Failed to merge $INP_NTUPLE_LIST and $INP_NANO_LIST" >> $(dirname $OUT_PATH)/failed.log
    echo $LOG_PATH
    rm $OUT_PATH
    exit 1
fi

rm $LOG_PATH
