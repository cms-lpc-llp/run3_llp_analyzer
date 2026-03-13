#!/bin/bash
# This script sets the environment variables for OpenSSL 1.1. Needed for EOS, etc, not automatically there for el9
#download SSL in some area ($HOME/openssl_libs) for me using "wget https://www.openssl.org/source/openssl-1.1.1w.tar.gz"
# and then unpack it with "tar -xvf openssl-1.1.1w.tar.gz"

export OPENSSL_HOME=$HOME/openssl_libs/openssl-1.1
export LD_LIBRARY_PATH=$OPENSSL_HOME/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
export PATH=$OPENSSL_HOME/bin:$PATH

# Remove CUDA stub path to avoid libnvidia-ml warnings.
CUDA_STUB_PATH="/cvmfs/cms.cern.ch/el8_amd64_gcc12/external/cuda/12.4.1-fc5cb0e72dba64b6abbf00089f3a044c/lib64/stubs"
if [ -n "${LD_LIBRARY_PATH}" ]; then
  new_ld=""
  old_ifs=$IFS
  IFS=:
  for p in $LD_LIBRARY_PATH; do
    if [ "$p" != "$CUDA_STUB_PATH" ]; then
      if [ -z "$new_ld" ]; then
        new_ld="$p"
      else
        new_ld="${new_ld}:$p"
      fi
    fi
  done
  IFS=$old_ifs
  export LD_LIBRARY_PATH="$new_ld"
fi
