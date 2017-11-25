#!/bin/bash

export CLHEP_PATH=$HOME/cern/clhep
#pushd .
#cd ~/cern
#source ./root/bin/thisroot.sh
#popd

mkdir -p build
cd build
cmake .. 
cmake --build . 
cd ..
