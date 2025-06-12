#!/bin/bash
# This script sets the environment variables for OpenSSL 1.1. Needed for EOS, etc, not automatically there for el9
#download SSL in some area ($HOME/openssl_libs) for me using "wget https://www.openssl.org/source/openssl-1.1.1w.tar.gz"
# and then unpack it with "tar -xvf openssl-1.1.1w.tar.gz"

export OPENSSL_HOME=$HOME/openssl_libs/openssl-1.1
export LD_LIBRARY_PATH=$OPENSSL_HOME/lib:$LD_LIBRARY_PATH
export PATH=$OPENSSL_HOME/bin:$PATH
