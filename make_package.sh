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
      doc/logo_ringmesh_doc.png
      doc/logo_ringmesh.png
      tests
      CMakeOptions.txt.sample
      CMakeLists.txt
      cmake"
       
tar vczf RINGMesh-$version.tar.gz $files
zip -r RINGMesh-$version.zip $files

