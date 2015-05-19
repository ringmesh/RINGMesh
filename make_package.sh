#!/bin/sh

if [ $# -eq 0 ]; then
    echo Enter a package version
    exit
else
    version=$1
fi

files="configure.bat
      configure.sh
      include
      src
      doc/Doxyfile
      doc/dox
      tests
      CMakeOptions.txt
      CMakeLists.txt
      cmake"
       
tar vczf RINGMesh-$version.tar.gz $files
zip -r RINGMesh-$version.zip $files

