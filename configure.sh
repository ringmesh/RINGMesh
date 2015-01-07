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

os="Linux64-gcc"
os_dynamic=$os-dynamic

#  Import plaform specific environment

. src/third_party/geogram/cmake/platforms/$os_dynamic/setvars.sh || exit 1

# Generate the Makefiles

echo
echo ================== GeoGram ====================

cat << EOF > src/third_party/geogram/CMakeOptions.txt
set(GEOGRAM_WITH_TETGEN TRUE)
set(GEOGRAM_WITH_MEDIT FALSE)
set(GEOGRAM_WITH_GRAPHICS FALSE)
EOF

for config in Release Debug
do
   platform=$os-$config
   echo
   echo ============= Creating makefiles for $platform ============
   build_dir=build/geogram/$platform

   mkdir -p $build_dir
   (cd $build_dir; $CMAKE -Wno-dev -DCMAKE_BUILD_TYPE:STRING=$config -DCMAKE_CXX_FLAGS:STRING="-fPIC" -DCMAKE_C_FLAGS:STRING="-fPIC" -DVORPALINE_PLATFORM:STRING=$os_dynamic ../../../src/third_party/geogram/; $CMAKE --build . -- -j4)
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

