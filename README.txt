Instructions for compiling RINGMesh

RINGMesh is tested under Linux (64 bits) and Windows (64 bits).
You will need CMake (version >= 2.8.11). There is no other dependancy (everything
 you need is shipped with RINGMesh). Follow the Linux or Windows instructions below.

Linux
=====

Configuring RINGMesh
------------------

Go to RINGMesh directory and then copy the file CMakeOptions.txt.sample into
CMakeOptions.txt: cp CMakeOptions.txt.sample CMakeOptions.txt.

If you need to use the graphic part of RINGMesh (visualization of RINGMesh::GeoModel
and RINGMesh::MacroMesh) set the variable RINGMESH_WITH_GRAPHICS to TRUE (uncomment the line).

Then to configure: ./configure.sh. This will first configure and compile Geogram 
and then it will generate dynamic libraries.

Compiling RINGMesh
------------------

- cd build/ringmesh/Linux64-gcc-Release
- make [-j4]

To build in debug, go to build/ringmesh/Linux64-gcc-Debug instead.

Eclipse-cdt project is provided (.project and .cproject). You can import RINGMesh into
Eclipse: File>Import...>General>Existing Projects into Workspace. Click on next, then
select the root directory (RINGMesh directory) then click on Finish. There are two
build configurations: release and debug. There is a target to clean and compile. There are
also targets to build/rebuild geogram.

Compiling the documentation
---------------------------

The documentation can be generated using [Doxygen](http://www.stack.nl/~dimitri/doxygen/):
- cd doc
- doxygen Doxyfile

Then you can open the doc/html/index.html file in your web browser. If you have qhelpgenerator
in your path, a .qch file is generated (ringmesh.qch). It can be imported in Qt Assisant: Edit>Preferences>
Documentation>add and select the .qch file.

Windows
=======

Configuring RINGMesh
------------------

Go to RINGMesh directory, open and save the CMakeOptions.txt.sample file as CMakeOptions.txt.
To use the graphic components of RINGMesh (visualization of RINGMesh::GeoModel
and RINGMesh::MacroMesh) set the variable RINGMESH_WITH_GRAPHICS to TRUE (uncomment the line).

To configure RINGMesh launch the configure.bat or by command line: call ./configure.bat.
By default, the project will be configured for Visual Studio 2010 in 64bits.
You can specify a different platform by adding it as an argument when executing
the configuration script in command line.
The supported platforms are :
- Visual Studio 10 2010 Win64
- Visual Studio 11 2012 Win64
- Visual Studio 12 2013 Win64


Compiling RINGMesh
------------------

You can either launch building in VisualStudio or calling cmake in command line
in the build directory created at the configuration step.
  cmake --build . --config Release


