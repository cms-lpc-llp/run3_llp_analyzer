#!/bin/bash


CHUNK_SIZE=1
TMP_PATH=/tmp/LLP_analyzer_tmpfiles
BIN=../bin/Runllp_MuonSystem_CA
N_JOBS=128
XRDCP=

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
        echo "  --xrdcp                 Xrdcp data local first. Only usable if CHUNK_SIZE=1"
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
    --xrdcp)
        XRDCP=1
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
echo "XRDCP          = $XRDCP"

if [ $CHUNK_SIZE -ne 1 ] && [ -n "$SYNC" ]; then
    echo "LOCAL and RSYNC cannot be used with CHUNK_SIZE > 1" 1>&2
    exit 1
fi

function get_label {
    path=$1
    YEAR=$(realpath $FILE | grep -o "Run[0-9]\{4\}.")
    YEAR=${YEAR:3:6}
    label=$(map_run_campaign Run$YEAR)
    if [ -n "$label" ]; then
        echo $label
        return
    fi

    if [[ $path == *"Summer22EE"* ]]; then
        echo "Summer22EE"
    elif [[ $path == *"Summer22"* ]]; then
        echo "Summer22"
    elif [[ $path == *"Summer23BPix"* ]]; then
        echo "Summer23BPix"
    elif [[ $path == *"Summer23"* ]]; then
        echo "Summer23"
    elif [[ $path == *"Summer24"* ]]; then
        echo "Summer24"
    else
        echo "Unknown campaign for $path" 1>&2
    fi
}

function map_run_campaign {
    run=$1 # Run20xxX
    if [ "$run" == "Run2022B" ] || [ "$run" == "Run2022C" ] || [ "$run" == "Run2022D" ]; then
        echo "Summer22"
    elif [ "$run" == "Run2022E" ] || [ "$run" == "Run2022F" ] || [ "$run" == "Run2022G" ]; then
        echo "Summer22EE"
    elif [ "$run" == "Run2023B" ] || [ "$run" == "Run2023C" ]; then
        echo "Summer23"
    elif [ "$run" == "Run2023D" ]; then
        echo "Summer23BPix"
    elif [[ "$run" == "Run2024"* ]]; then
        echo "Summer24"
    fi
}

function prepare_chunks {
    
    FILE=$1
    LOCAL_TMP_PATH=$2/$(basename ${FILE%.txt})
    CHUNK_SIZE=$3

    mkdir -p $LOCAL_TMP_PATH

    i=0
    chunk=()
    while IFS= read -r line
    do
        if [[ $line != *.root ]]; then
            continue
        fi
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
    YEAR=$(realpath $FILE | grep -o "Run[0-9]\{4\}.")
    YEAR=${YEAR:3:6}

    LABEL=$(get_label $FILE)
    

    echo -e "YEAR=$YEAR\nIS_DATA=$IS_DATA\nLABEL=$LABEL\n" > $LOCAL_TMP_PATH/config.env
}


function launch {
    LIST_FILE=$1
    OUT_PATH=$2
    BIN=$3
    XRDCP=$4
    # read config from the same path, contains IS_DATA, LABEL, YEAR
    source $(dirname $LIST_FILE)/config.env
    FILE_OUT_PATH=$OUT_PATH/$(basename $(dirname $LIST_FILE))/$(basename ${LIST_FILE%.txt}).root

    if [ -f $FILE_OUT_PATH ]; then
        return
    fi

    if [ -n "$XRDCP" ]; then
        f=$(cat $LIST_FILE)
        LOCAL_ROOT_PATH=${LIST_FILE%.txt}.root
        xrdcp -f $f $LOCAL_ROOT_PATH > ${LIST_FILE%.txt}.sync.log 2>&1
        echo $LOCAL_ROOT_PATH > $LIST_FILE
    fi


    mkdir -p $(dirname $FILE_OUT_PATH)
    $BIN $LIST_FILE -d=$IS_DATA -l=$LABEL -f=$FILE_OUT_PATH > ${LIST_FILE%.txt}.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "Failed to run $LIST_FILE: log at ${LIST_FILE%.txt}.log" 1>&2
        rm $FILE_OUT_PATH
    else
        rm $LIST_FILE
        rm ${LIST_FILE%.txt}.log
        if [ -n "$LOCAL_ROOT_PATH" ]; then
            rm $LOCAL_ROOT_PATH
            rm ${LIST_FILE%.txt}.sync.log
        fi
    fi
}

export -f launch

for INP_LIST in $INP_LIST_FILES; do
    prepare_chunks $INP_LIST $TMP_PATH $CHUNK_SIZE
done

echo $TMP_PATH/*/*.txt | tr ' ' '\n' | parallel --lb -j $N_JOBS --bar --eta launch {} $OUT_PATH $BIN $XRDCP
