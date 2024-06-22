#!/bin/bash
SRC=$1
DST_DIR=$2

if [[ "$SRC" == "/eos/uscms/"* ]]; then
    DST=$(echo $SRC | sed "s|/eos/uscms|$DST_DIR|")
    mkdir -p $(dirname $DST)
    if [ -f $DST ]; then
        exit 0
    fi

    rsync -a lpc:$SRC $DST 
fi

DST=$(echo $SRC | sed "s|root://cmsxrootd.fnal.gov/|$DST_DIR|")
echo $DST
if [ -f $DST ]; then
    exit 0
fi

mkdir -p $(dirname $DST)
echo xrdcp $SRC $DST --retry 3 --rm-bad-cksum
xrdcp $SRC $DST --retry 3 --rm-bad-cksum > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Failed to copy $SRC" >> log.err
    exit 1
fi
