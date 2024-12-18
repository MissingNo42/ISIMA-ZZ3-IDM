#!/usr/bin/env sh

git clone https://gitlab.cern.ch/CLHEP/CLHEP.git

mkdir build_clhep
cd build_clhep
cmake ../CLHEP
cmake --build . --config RelWithDebInfo
ctest
cmake --build . --target install
