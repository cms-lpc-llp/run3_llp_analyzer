#!/bin/sh

outpath=$1
shift
njobs=$1
shift

function worker() {
    inp_path=$1
    out_path=$2
    out_name=$(echo $inp_path | sed -r "s/.+\/(.+)\/(.+)\/(.+)\/(.+)\..+/\1_\2\3_\4/")
    if [ -f $out_path/$out_name.h5 ]; then
        return
    fi

    if [[ $inp_path == "/eos/uscms"* ]]; then
        inp_path="root://cmseos.fnal.gov/"${inp_path#/eos/uscms}    
    fi

    if [[ $inp_path == "root://cmseos.fnal.gov/"* ]]; then
        inp_path="/storage/cms"${inp_path#root://cmseos.fnal.gov/}    
    fi


    if [[ $inp_path == "root://cmsxrootd.fnal.gov/"* ]]; then
        inp_path="/storage/cms"${inp_path#root://cmsxrootd.fnal.gov/}
    fi
    #echo $inp_path
    #if [ ! -f $inp_path ]; then
    #    echo $inp_path not found
    #fi
    python cache_one.py -i $inp_path -o $out_path/.$out_name.h5

    succeed=$?
    if [ $succeed -eq 0 ]; then
        mv $out_path/.$out_name.h5 $out_path/$out_name.h5
    fi
}

export -f worker

for list in $@; do
    mkdir -p $outpath/$(basename $list)
    cat $list | parallel --progress -j $njobs worker {} $outpath/$(basename $list)
done
