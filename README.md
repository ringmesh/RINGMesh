[![Travis build](https://ci.appveyor.com/api/projects/status/nlso0s96wcuge2vn/branch/master?svg=true)](https://ci.appveyor.com/project/ringmesh/ringmesh/branch/master)
[![Appveyor build](https://travis-ci.org/ringmesh/RINGMesh.svg?branch=master)](https://travis-ci.org/ringmesh/RINGMesh)
[![Sonar quality](https://sonarcloud.io/api/badges/gate?key=ringmesh)](https://sonarcloud.io/dashboard/index/ringmesh)
[![Coverage](https://sonarcloud.io/api/badges/measure?key=ringmesh&metric=coverage)](https://sonarcloud.io/dashboard/index/ringmesh)

# Configuration and compilation

RINGMesh is tested under Linux (64 bits) and Windows (64 bits).
You will need CMake (version >= 3.1). There is no other dependency (everything
you need is shipped with RINGMesh). Follow the Linux, Mac OS or Windows instructions below.

To clone RINGMesh you need Git (https://git-scm.com/).
Make sure that Git binary directory is in your computer path (environment variable).
Under Windows, after installing Git you should have in your path environment variable:
C:\Program Files\Git\cmd.
Warning: TortoiseGit (https://tortoisegit.org/) does not install Git.

You can clone RINGMesh using:

```
git clone https://github.com/ringmesh/RINGMesh/
```

## Linux


### Configuring RINGMesh

Execute cmake command in a RINGMesh/build directory.

```bash
mkdir build
cd build
```
To configure using default options:
```bash
cmake ..
```
To define the options, use the cmake interface

```bash
cmake-gui .. or ccmake ..
```
 Note: Configuration options can also be edited by creating a `UserConfig.cmake` file from a `DefaultConfig.cmake` file
 in the cmake folder.

### Compiling RINGMesh

To compile RINGMesh, go to RINGMesh root directory and:

```bash
cd build/Release
make
```
To build in debug, go to build/Debug instead.

Note: RINGMesh uses C++ 11 features. You need gcc/g++ version higher or equal to 4.8 to compile it.

### Generating the documentation

If Doxygen is installed on your computer, a target ```doc-devkit``` is built during RINGMesh configuration. You can generate the documentation using the following command:

```bash
make doc-devkit
```

### Troubleshooting

If you get this error during geogram gfx compilation (occured for Ubuntu 17.04):
```
No rule to make target '/usr/lib/x86_64-linux-gnu/libGL.so'
```
you can try under root:
```bash
rm /usr/lib/x86_64-linux-gnu/libGL.so
ln -s /usr/lib/libGL.so.1 /usr/lib/x86_64-linux-gnu/libGL.so
```
### Additionnal information

#### Eclipse-cdt project
[Eclipse-cdt](http://www.eclipse.org/cdt/)
project is provided (.project and .cproject). You can import RINGMesh into
Eclipse: File>Import...>General>Existing Projects into Workspace.

## Windows

### Configuring RINGMesh

 * Launch CMake GUI interface. 
 * Indicate where is the source code as the path to RINGMesh root and where to put the binaries as 
 `RINGMesh_root/build`.
 * Set the Configuration options using the CMake GUI interface.
 * Launch the `configure`and `generate` option.

 Note: Configuration options can also be edited by creating a `UserConfig.cmake` file from a `DefaultConfig.cmake` file
 in the cmake folder.

RINGMesh compils with the following visual studio version:

* Visual Studio 12 2013 Win64
* Visual Studio 14 2015 Win64
* Visual Studio 15 2017 Win64

Note: RINGMesh uses C++11 features. Make sure that you have installed C++ package for VisualStudio through the 
VisualStudio installer.


### Compiling RINGMesh

RINGMesh need several third parties that are automatically compiled and installed by compiling 
the project `SUPERBUILD.sln`.

 * Open the project `SUPERBUILD.sln` with the visual studio version that you choose during the configuration step.
 * Build the solution.
 
All the RINGMesh third parties have now been compiled, installed and the `RINGMesh.sln` project have been created.

 * Open the project `RINGMesh.sln`.
 * Build the solution. 

The available compilation modes are:

* Release
* Debug
* RelWithDebInfo 

Note: The compilation mode `RelWithDebInfo` is mandatory to run a Gocad plugin in Debug mode with Gocad
  in Release, there are issues between libraries in Debug linked to a Gocad plugin.

### Compiling the documentation

* Check the BUILD_DOCUMENTATION option during the configuration of RINGMesh with cmake
* Open the solution which is in build/RINGMesh.sln in VisualStudio
* Build the doc-devkit

See the documentation section for more details.

## Mac OS

### Configuring RINGMesh

#### Using clang (without Xcode)
As in Linux.

#### Using Xcode IDE
As in Windows but with the Xcode generator
(use "-G Xcode" if you use cmake in command lines).

### Compiling RINGMesh
You need to install the Mac OS "Command Line Developer Tools".

Note: you need gcc/g++ version higher or equal to 4.2 to compile RINGMesh.
In Mac OS, clang is used.

#### Using clang (without Xcode)
As in Linux except for the packages.

#### Using Xcode IDE
You need to install Xcode IDE.
Open the build/ringmesh/RINGMesh.xcodeproj with Xcode IDE,
and then compile (as in Windows with VisualStudio).
Or use these command lines:
```
cd build/ringmesh
xcodebuild -project RINGMesh.xcodeproj -alltargets -configuration Release
```
To build in Debug, replace "Release" by "Debug" after "-configuration".
