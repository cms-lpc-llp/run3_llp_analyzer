#!/bin/bash


CHUNK_SIZE=1
TMP_PATH=/tmp/LLP_analyzer_tmpfiles
BIN=../bin/Runllp_MuonSystem_CA
N_JOBS=128
SYNC=

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
        echo "Usage: $0 [options] -- [input list files]"
        echo "Options:"
        echo "  -o, --output <path>     Output path"
        echo "  -t, --temp <path>       Temporary path, default /tmp/LLP_analyzer_tmpfiles"
        echo "  -j, --jobs <number>     Number of parallel jobs, default 128"
        echo "  -c, --chunk <number>    Number of files per job, default 1"
        echo "  -b, --bin <path>        Path to the binary, default ../bin/Runllp_MuonSystem_CA"
        echo "  --rsync <path>          Rsync path to the data, default none. Only usable if CHUNK_SIZE=1"
        echo "  --local                 Use local path instead of xrootd. Use this if you are on LPC. Only usable if CHUNK_SIZE=1"
        exit 0
        ;;
    -o|--output)
        OUT_PATH="$2"
        shift
        shift
        ;;
    -t|--temp)
        TMP_PATH="$2"
        shift
        shift
        ;;
    -j|--jobs)
        N_JOBS="$2"
        shift
        shift
        ;;
    -c|--chunk)
        CHUNK_SIZE="$2"
        shift
        shift
        ;;
    -b|--bin)
        BIN="$2"
        shift
        shift
        ;;
    --rsync)
        SYNC="$2"
        shift
        shift
        ;;
    --local)
        SYNC="_local_"
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

echo "INP_LIST_FILES  = $INP_LIST_FILES"
echo "OUT_PATH       = $OUT_PATH"
echo "CHUNK_SIZE     = $CHUNK_SIZE"
echo "TMP_PATH       = $TMP_PATH"
echo "N_JOBS         = $N_JOBS"
echo "BIN            = $BIN"
echo "SYNC           = $SYNC"

if [ $CHUNK_SIZE -ne 1 ] && [ -n "$SYNC" ]; then
    echo "LOCAL and RSYNC cannot be used with CHUNK_SIZE > 1" 1>&2
    exit 1
fi


function prepare_chunks {
    
    FILE=$1
    LOCAL_TMP_PATH=$2/$(basename ${FILE%.list})
    CHUNK_SIZE=$3
    LABEL=test

    mkdir -p $LOCAL_TMP_PATH

    i=0
    chunk=()
    while IFS= read -r line
    do
        chunk+=("$line")
        if (( ${#chunk[@]} == CHUNK_SIZE )); then
            i=$((i+1))
            printf "%s\n" "${chunk[@]}" > $LOCAL_TMP_PATH/$i.txt
            chunk=()
        fi
    done < "$FILE"
    if (( ${#chunk[@]} > 0 )); then
        i=$((i+1))
        printf "$i %s\n" "${chunk[@]}" > $LOCAL_TMP_PATH/$i.txt
    fi

    IS_DATA=$(realpath $FILE | grep -q '\(202[0-9][A-Z]\|Data\)' && echo "yes" || echo "no")
    YEAR=$(realpath $FILE | grep -o "Data[0-9]\{4\}")
    YEAR=${YEAR:4:4}

    echo -e "YEAR=$YEAR\nIS_DATA=$IS_DATA\nLABEL=$LABEL\n" > $LOCAL_TMP_PATH/config.env
}


function launch {
    LIST_FILE=$1
    OUT_PATH=$2
    BIN=$3
    SYNC=$4
    # read config from the same path, contains IS_DATA, LABEL, YEAR
    source $(dirname $LIST_FILE)/config.env
    FILE_OUT_PATH=$OUT_PATH/$(basename $(dirname $LIST_FILE))/$(basename ${LIST_FILE%.list}).root

    if [ -f $FILE_OUT_PATH ]; then
        return
    fi

    if [ -n "$SYNC" ]; then
        prefix="root://cmsxrootd.fnal.gov/"
        f=$(cat $LIST_FILE)
        if [ "$SYNC" == "_local_" ]; then
            LOCAL_ROOT_PATH=/eos/uscms/${f#$prefix}
        else
            LOCAL_ROOT_PATH=${LIST_FILE%.list}.root
            if [ ! -f $LOCAL_ROOT_PATH ]; then
                rsync -a $SYNC/${f#$prefix} $LOCAL_ROOT_PATH > /dev/null 2> >(grep -i -v 'Warning: Permanently')
            fi
        fi
        echo $LOCAL_ROOT_PATH > $LIST_FILE
    fi


    mkdir -p $(dirname $FILE_OUT_PATH)
    $BIN $LIST_FILE -d=$IS_DATA -l=$LABEL -f=$FILE_OUT_PATH > ${LIST_FILE%.list}.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "Failed to run $LIST_FILE: log at ${LIST_FILE%.list}.log" 1>&2
        rm $FILE_OUT_PATH
    else
        rm $LIST_FILE
        rm ${LIST_FILE%.list}.log
        if [ -n "$SYNC" ] && [ "$SYNC" != "_local_" ]; then
            rm $LOCAL_ROOT_PATH
        fi
    fi
}

export -f launch

for INP_LIST in $INP_LIST_FILES; do
    prepare_chunks $INP_LIST $TMP_PATH $CHUNK_SIZE
done

echo $TMP_PATH/*/*.txt | tr ' ' '\n' | parallel -j $N_JOBS --bar --eta launch {} $OUT_PATH $BIN $SYNC
