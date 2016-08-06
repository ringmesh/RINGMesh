#!/bin/bash

# Get cmake binaries
wget --no-check-certificate https://cmake.org/files/v3.5/cmake-3.5.2-Linux-x86_64.tar.gz
tar -xzf cmake-3.5.2-Linux-x86_64.tar.gz

# Create build tree
mkdir build
cd build

# Configure & build & test
../cmake-3.5.2-Linux-x86_64/bin/cmake ..
cd ringmesh/Debug
make
make test
