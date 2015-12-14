#!/bin/sh

if [ $# -eq 0 ]; then
    echo Enter a package version
    exit
else
    version=$1
fi


cd `dirname $0`
RINGMesh_HOME=$(pwd)

files="configure.bat
      configure.sh
      include
      src
      doc/Doxyfile
      doc/dox
      doc/logo_ringmesh_doc.png
      doc/logo_ringmesh.png
      doc/ringmesh.png
      doc/bme_doc.png
      tests
      CMakeOptions.txt.sample
      CMakeLists.txt
      cmake"

main_directory=RINGMesh
if [ -d ${main_directory} ]
then
    rm -rf ${main_directory}
fi
mkdir ${main_directory}
cp -R --parents ${files} ${main_directory}
       
tar vczf RINGMesh-${version}.tar.gz ${main_directory}
zip -r RINGMesh-${version}.zip ${main_directory}

rm -rf ${main_directory}

doc_target=doc_online
cd doc
rm -rf $doc_target
doxygen Doxyfile_online
cd $doc_target/html
tar vczf ${RINGMesh_HOME}/doc.tar.gz ./*


