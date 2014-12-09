#!/bin/sh
if [ $# -eq 0 ]; then
    CMAKE=cmake
else
    CMAKE=$1
fi

echo ============= Checking for cmake ============
if ($CMAKE --version)
then
   echo "found CMake"
else
   echo "cmake not found, please install it (see http://www.cmake.org/)" 
   exit
fi

os="$1"
if [ -z "$os" ]; then
    os=`uname -a`
    case "$os" in
        Linux*x86_64*)
            os=Linux64-gcc
            ;;
        Linux*amd64*)
            os=Linux64-gcc
            ;;
        Linux*i586*|Linux*i686*)
            os=Linux32-gcc
            ;;
        *)
            echo "Error: OS not supported: $os"
            exit 1
            ;;
    esac
fi

#  Import plaform specific environment

. src/third_party/geogram/cmake/platforms/$os/setvars.sh || exit 1

# Generate the Makefiles

echo
echo ================== GeoGram ====================

for config in Release Debug
do
   platform=$os-$config
   echo
   echo ============= Creating makefiles for $platform ============
   build_dir=build/geogram/$platform

   mkdir -p $build_dir
   (cd $build_dir; $CMAKE -Wno-dev -DCMAKE_BUILD_TYPE:STRING=$config -DCMAKE_CXX_FLAGS:STRING="-fPIC" -DGEOGRAM_WITH_TETGEN:BOOL=TRUE -DCMAKE_C_FLAGS:STRING="-fPIC" -DVORPALINE_PLATFORM:STRING=$os ../../../src/third_party/geogram/; $CMAKE --build .)
done
echo
echo ============== GeoGram build configured ==================
echo


echo
echo ================== GRGMesh ====================

for config in Release Debug
do
   platform=$os-$config
   echo
   echo ============= Creating makefiles for $platform ============
   build_dir=build/grgmesh/$platform

   mkdir -p $build_dir
   (cd $build_dir; $CMAKE -Wno-dev -DCMAKE_BUILD_TYPE:STRING=$config -DGEOGRAM_PLATFORM:STRING=$os ../../../)
done
echo
echo ============== GRGMesh build configured ==================

cat << EOF
to build:
  go to build/$os-Release or build/$os-Debug
  and 'make'
Note: local configuration can be specified in CMakeOptions.txt
(see CMakeOptions.txt.sample for an example)
You'll need to re-run configure.sh if you create or modify CMakeOptions.txt
EOF

