#!/bin/sh

cd /afs/cern.ch/work/f/fernanpe/CMSSW_10_6_32/src/run3_llp_analyzer/DNN
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

python reader.py --inputFile $1
