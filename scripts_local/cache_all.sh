#!/bin/sh

outpath=$1
shift
njobs=$1
shift

function worker() {
    inp_path=$1
    out_path=$2
    out_name=$(echo $inp_path | sed -r "s/.+\/(.+)\/(.+)\/(.+)\/(.+)\..+/\1_\2\3_\4/")
    if [ -f $out_path/$out_name ]; then
        return
    fi

    if [ $inp_path == "/eos/uscms"* ]; then
        inp_path="root://cmseos.fnal.gov/"${inp_path#/eos/uscms}    
    fi

    python cache_one.py -i $inp_path -o $out_path/$out_name.tmp

    succeed=$?
    if [ $succeed -eq 0 ]; then
        mv $out_path/$out_name.tmp $out_path/$out_name
    fi
}

export -f worker

for list in $@; do
    mkdir -p $outpath/$(basename $list)
    cat $list | parallel --progress -j $njobs worker {} $outpath/$(basename $list)
done