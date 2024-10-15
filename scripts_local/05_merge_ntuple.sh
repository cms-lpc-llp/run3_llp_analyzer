#!/bin/bash

CACHE=./cache
N_JOBS=128
BIN=../MergeNtuples
BSZ_NTUPLE=50
TMP_PATH=/tmp/MergeNtuples
BSZ_NANO=5
NANOS=""
NTUPLES=""
OUT_DIR=./output

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
        echo "Usage: $0 [options] <input_list_files>"
        echo "Options:"
        echo "  -c, --cache <path>              Path for storing the output caches, defaults to ./cache"
        echo "  -j, --jobs <n>                  Number of parallel jobs, defaults to 128"
        echo "  -b, --bin <path>                Path to the binary, defaults to ../CacheNtuples"
        echo "  -B1, --batch_ntuple <n>         Number of ntuple files to merge in a single job, defaults to 1"
        echo "  -B2, --batch_nano <n>           Number of nano files to merge in a single job, defaults to 1"
        echo "  -t, --tmp_path <path>           Path for storing the temporary files, defaults to /tmp/MergeNtuples"
        echo "  -nt, --ntuple <paths>           Paths to the ntuple files"
        echo "  -na, --nano <paths>             Paths to the nano files"
        echo "  -o, --out_dir <path>            Path for storing the output files, defaults to ./output"
        echo "  -l, --local <path>              Replace path prefix by this, useful if everything presents locally"
        exit 0
        ;;
    -B1|--batch_ntuple)
        BSZ_NTUPLE="$2"
        shift
        shift
        ;;
    -B2|--batch_nano)
        BSZ_NANO="$2"
        shift
        shift
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
    -l|--local)
        LOCAL="$2"
        shift
        shift
        ;;
    -o|--out-dir)
        OUT_DIR="$2"
        shift
        shift
        ;;
    -na|--nano)
        shift
        while [[ $# -gt 0 ]]; do
            case $1 in
            -*|--*)
                break
                ;;
            *)
                NANOS="$NANOS $1"
                shift
                ;;
            esac
        done
        ;;
    -nt|--ntuple)
        shift
        while [[ $# -gt 0 ]]; do
            case $1 in
            -*|--*)
                break
                ;;
            *)
                NTUPLES="$NTUPLES $1"
                shift
                ;;
            esac
        done
        ;;
    -t|--tmp_path)
        TMP_PATH="$2"
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

echo "CACHE = $CACHE"
echo "N_JOBS = $N_JOBS"
echo "BIN = $BIN"
echo "TMP_PATH = $TMP_PATH"
echo "BSZ_NTUPLE = $BSZ_NTUPLE"
echo "BSZ_NANO = $BSZ_NANO"
echo ================================================================================
echo "NTUPLES = $NTUPLES"
echo ================================================================================
echo "NANOS = $NANOS"
echo ================================================================================

# if [ ! -f $BIN ]; then
#   echo "Binary $BIN not found"
#   exit 1
# fi


function prepare_chunks {
    
    FILE=$(realpath $1)
    TMP_PATH=$(realpath $2)
    CHUNK_SIZE=$3
    TYPE=$4
    LOCAL=$5

    YEAR_ERA=$(echo $FILE | grep -oP "(?<=Run)[0-9]{4}[A-Z]")
    VERSION=$(echo $FILE | grep -oP '(?<=(\-|_)v)(?<!v[0-9]\-v)[0-9]')

    LOCAL_TMP_PATH=$TMP_PATH/$YEAR_ERA/$TYPE
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
        printf "%s\n" "${chunk[@]}" > $LOCAL_TMP_PATH/$i.txt
    fi

    if [ ! -z $LOCAL ]; then
        sed -i "s|/eos/uscms|$LOCAL|" $LOCAL_TMP_PATH/*.txt
        sed -i "s|root://cmsxrootd.fnal.gov/|$LOCAL|" $LOCAL_TMP_PATH/*.txt
    fi
    
}

function argsgen {
    TMP_PATH=$(realpath $1)
    for dir in $(ls $TMP_PATH); do
        if [ ! -d $TMP_PATH/$dir/ntuple ] || [ ! -d $TMP_PATH/$dir/nano ]; then
            continue
        fi
        for ntuple_f in $(ls $TMP_PATH/$dir/ntuple); do
            for nano_f in $(ls $TMP_PATH/$dir/nano); do
                echo "$TMP_PATH/$dir/ntuple/$ntuple_f $TMP_PATH/$dir/nano/$nano_f"
            done
        done
    done
}

function launch {
    BIN=$1
    LIST_NTUPLE=$2
    LIST_NANO=$3
    CACHE_PATH=$4
    OUT_DIR=$5


    YEAR_ERA_VERSION=$(basename $(dirname $(dirname $LIST_NTUPLE)))
    FNAME="$(basename $LIST_NTUPLE .txt)-$(basename $LIST_NANO .txt)"
    
    mkdir -p $OUT_DIR/$YEAR_ERA_VERSION
    OUT_PATH="$OUT_DIR/$YEAR_ERA_VERSION/$FNAME.root"

    if [ -f $OUT_PATH ] || [ -f $OUT_PATH.nomatch ]; then
        return
    fi

    $BIN $LIST_NTUPLE $LIST_NANO $OUT_PATH $CACHE_PATH > $OUT_DIR/$YEAR_ERA_VERSION/$FNAME.log 2>&1

    if [ $? -ne 0 ]; then
        echo "Failed to merge $LIST_NTUPLE and $LIST_NANO" 1>&2
        echo "Fauled to merge $LIST_NTUPLE and $LIST_NANO" >> $OUT_DIR/failed.log
        rm $OUT_PATH
        exit 1
    fi


    if [ ! -f $OUT_PATH ]; then
        touch $OUT_PATH.nomatch
    fi

    rm $OUT_DIR/$YEAR_ERA_VERSION/$FNAME.log
}


for NANO_LIST in $NANOS; do
    prepare_chunks $NANO_LIST $TMP_PATH $BSZ_NANO nano $LOCAL
done

for NTUPLE_LIST in $NTUPLES; do
    prepare_chunks $NTUPLE_LIST $TMP_PATH $BSZ_NTUPLE ntuple $LOCAL
done

export -f launch

argsgen $TMP_PATH | parallel -j $N_JOBS --colsep ' ' --eta --bar launch $BIN {1} {2} $CACHE $OUT_DIR