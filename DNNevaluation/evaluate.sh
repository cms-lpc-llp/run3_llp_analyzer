#!/bin/sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh

cd /afs/cern.ch/work/f/fernanpe/CMSSW_10_6_32/src/run3_llp_analyzer/DNNevaluation_240310
python EvaluateDNN.py --in_file $1
