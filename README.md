How to compile RINGMesh             {#ringmesh_compiling}
=======================

[![Appveyor build](https://ci.appveyor.com/api/projects/status/nlso0s96wcuge2vn/branch/master?svg=true)](https://ci.appveyor.com/project/ringmesh/ringmesh/branch/master)
[![Travis build](https://travis-ci.org/ringmesh/RINGMesh.svg?branch=master)](https://travis-ci.org/ringmesh/RINGMesh)
[![Sonar quality](https://sonarcloud.io/api/badges/gate?key=ringmesh)](https://sonarcloud.io/dashboard/index/ringmesh)
[![Coverage](https://sonarcloud.io/api/badges/measure?key=ringmesh&metric=coverage)](https://sonarcloud.io/dashboard/index/ringmesh)

RINGMesh is tested under Linux (64 bits) and Windows (64 bits).
You will need CMake (version >= 3.1). There is no other dependency (everything
you need is shipped with RINGMesh). Follow the Linux, Mac OS or Windows instructions below.

To clone RINGMesh you need Git (https://git-scm.com/).
Make sure that Git binary directory is in your computer path (environment variable).
Under Windows, after installing Git you should have in your path environment variable:
C:\Program Files\Git\cmd.
Warning: TortoiseGit (https://tortoisegit.org/) does not install Git.

Linux
===========================

Configuring RINGMesh
--------------------

Execute cmake command in a RINGMesh/build directory.

* mkdir build
* cd build

To configure using default options:

* cmake ..

To define the options, use the cmake interface:

* cmake-gui .. or ccmake ..


Compiling RINGMesh
------------------

To compile you need the following packages (on Debian-based linux):
* build-essential
* libx11-dev
* libxrandr-dev
* libxinerama-dev
* libxcursor-dev
* freeglut3-dev
* libxi-dev

Note: you need gcc/g++ version higher or equal to 4.8 to compile RINGMesh.

Then, to compile RINGMesh, go to RINGMesh root directory and:

* cd build/Release
* make [-j4]

To build in debug, go to build/Debug instead.

Note: if you get this error during geogram gfx compilation
"No rule to make target '/usr/lib/x86_64-linux-gnu/libGL.so'"
(occured for Ubuntu 17.04), do:
* sudo rm /usr/lib/x86_64-linux-gnu/libGL.so
* sudo ln -s /usr/lib/libGL.so.1 /usr/lib/x86_64-linux-gnu/libGL.so

Eclipse-cdt project is provided (.project and .cproject). You can import RINGMesh into
Eclipse: File>Import...>General>Existing Projects into Workspace. Click on next, then
select the root directory (RINGMesh directory) then click on Finish. There are two
build configurations: release and debug. There is a target to clean and compile. There are
also targets to build/rebuild geogram.

Compiling the documentation
---------------------------

* Check the BUILD_DOCUMENTATION option when using cmake
 * cd build
 * ccmake ..
 * set BUILD_DOCUMENTATION option to ON
 * configure and generate
* cd build/Release
* make doc-devkit OR make doc-devkit-lite

You can also build the documentation through eclipse (see available targets).
See the documentation section for more details.

Windows
=======

Configuring RINGMesh
--------------------

Launch CMake GUI, indicate where is the source code as the path to RINGMesh root and
where to put the binaries as this_root/build.
Configuration options can be set in using the interface.

RINGMesh has previously been compiled with:

* Visual Studio 10 2010 Win64
* Visual Studio 11 2012 Win64
* Visual Studio 12 2013 Win64
* Visual Studio 14 2015 Win64
* Visual Studio 15 2017 Win64

Compiling RINGMesh
------------------

Make sure that you have installed C++ package for VisualStudio through the VisualStudio installer.
You can either launch building in VisualStudio or calling cmake in command line
in the build directory created at the configuration step:

* cmake --build . --config Release
* cmake --build . --config Debug
* cmake --build . --config RelWithDebInfo

The available compilation modes are:

* Release
* Debug
* RelWithDebInfo (mandatory to debug a Gocad plugin in Debug mode with a Gocad
  in Release, there are issues between libraries in Debug linked to a Gocad plugin)

Compiling the documentation
---------------------------

* Check the BUILD_DOCUMENTATION option when using cmake
* Open the solution which is in build/RINGmesh.sln in VisualStudio
* Build the doc-devkit or the doc-devkit-lite project

See the documentation section for more details.

Mac OS
======

Configuring RINGMesh
--------------------
### Using clang (without Xcode)
As in Linux.

### Using Xcode IDE
As in Windows but with the Xcode generator
(use "-G Xcode" if you use cmake in command lines).

Compiling RINGMesh
------------------
You need to install the Mac OS "Command Line Developer Tools".

Note: you need gcc/g++ version higher or equal to 4.2 to compile RINGMesh.
In Mac OS, clang is used.

### Using clang (without Xcode)
As in Linux except for the packages.

### Using Xcode IDE
You need to install Xcode IDE.
Open the build/RINGMesh.xcodeproj with Xcode IDE,
and then compile (as in Windows with VisualStudio).
Or use these command lines:
* cd build
* xcodebuild -project RINGMesh.xcodeproj -alltargets -configuration Release

To build in Debug, replace "Release" by "Debug" after "-configuration".

Compiling the documentation
---------------------------

See the documentation section for more details.

### Using clang (without Xcode)
As in Linux.

### Using Xcode IDE
As in Windows with VisualStudio but with Xcode.

About documentation
===================

The documentation can be generated using [Doxygen](http://www.stack.nl/~dimitri/doxygen/).
Two targets are available:
* doc-devkit include full documentation of RINGMesh and Geogram
* doc-devkit-lite include only the RINGMesh documentation

Then you can go in doc/devkit[-lite]/html and open the index.html with your web browser.
A ringmesh.qch file is generated in doc/devkit[-lite]/html (to load in Qt Assistant).
