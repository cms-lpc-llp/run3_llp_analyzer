#!/bin/bash

version=3.4.2
TMP_PATH=$(mktemp -d)
REPO_PATH=$(pwd)

if [[ ! -n $(head -n1 README.md | grep run3_llp_analyzer) ]]; then
    echo "This script must be run from the run3_llp_analyzer directory" 1>&2
    exit 1
fi

if [ -d ./fastjet-install ]; then
    echo "FastJet is already installed at fastjet-install" 1>&2
    echo "If you wish to reinstall, please remove the fastjet-install directory" 1>&2
    exit 1
fi

echo "Downloading FastJet"
wget https://fastjet.fr/repo/fastjet-$version.tar.gz -O $TMP_PATH/fastjet-$version.tar.gz
echo b3d33155b55ce43f420cd6d99b525acf7bdc2593a7bb7ea898a9ddb3d8ca38e3 $TMP_PATH/fastjet-$version.tar.gz | sha256sum -c

if [ $? -ne 0 ]; then
    echo "FastJet download failed." 1>&2
    echo "Your file is kept at $TMP_PATH/fastjet-$version.tar.gz" 1>&2
    exit 1
fi

cd ..
tar -xzf $TMP_PATH/fastjet-$version.tar.gz

echo "Installing FastJet"
cd fastjet-$version

./configure --prefix="$REPO_PATH/fastjet-install"

if [ $? -ne 0 ]; then
    echo "FastJet configure failed"
    rm -rf $TMP_PATH
    rm -rf $REPO_PATH/../fastjet-$version
    exit 1
fi

make "$@" && make check && make install

cd -
rm -rf $TMP_PATH
echo FastJet installed at "$REPO_PATH/fastjet-install"
