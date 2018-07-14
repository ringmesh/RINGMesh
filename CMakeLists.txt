# Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
# Applications (ASGA). All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of ASGA nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#     http://www.ring-team.org
#
#     RING Project
#     Ecole Nationale Superieure de Geologie - GeoRessources
#     2 Rue du Doyen Marcel Roubault - TSA 70605
#     54518 VANDOEUVRE-LES-NANCY
#     FRANCE
#------------------------------------------------------------------------------------------------
# Root CMakeList of RINGMesh project
#------------------------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.1)

# Configuration modes - Debug and Release by default + RelWithDebInfo for Windows.
# RelWithDebInfo is mandatory to debug some plugins in Windows in debug mode (e.g., a Gocad plugin) to
# avoid issues in particular in the STL (e.g., bad transfer of std::string from Gocad to RINGMesh).
set(CMAKE_CONFIGURATION_TYPES Debug Release RelWithDebInfo CACHE TYPE INTERNAL FORCE)

#------------------------------------------------------------------------------------------------
# Options of the RINGMesh project
# They can be set
# 1 - using CMake GUI,
# 2 - in command line: cmake -D<var>:<type>=<value> MyProjectFolder
# Project options
# Optional components of RINGMesh - creation of specific targets
option(RINGMESH_WITH_GRAPHICS "Compile viewer" OFF)
option(RINGMESH_WITH_UTILITIES "Compile utility executables" OFF)
option(RINGMESH_WITH_TUTORIALS "Compile API trainings and tuturials" OFF)
option(RINGMESH_WITH_TESTS "Compile test projects" OFF)
option(RINGMESH_WITH_GUI "Compile GUI bindings" OFF)
option(RINGMESH_WITH_METIS "Unse RINGmesh with METIS for GeoModel partitionning" OFF)

# Optional custom steps 
option(BUILD_GEOGRAM_WITHOUT_EXE "Compile geogram without tests and executables" ON)
option(USE_SUPERBUILD "Whether or not a superbuild should be invoked" ON)
option(EXPORT_COMPILE_COMMANDS "Whether or not to build the compile commands in a json file" ON)
mark_as_advanced(BUILD_GEOGRAM_WITHOUT_EXE USE_SUPERBUILD RINGMESH_WITH_GUI RINGMESH_WITH_METIS)

# Select Available Mesher 
# Each selected mesher will increase the RINGMesh functionalities.
# By default RINGMesh is compiled with Tetgen 
# Check Tetgen licences to ensure you have the right to use it freely
option(RINGMESH_WITH_TETGEN "Use Tetgen tetrahedral mesher" ON)
# If a path is specified RINGMesh will be compiled to work with - empty by default
set(MG_TETRA "" CACHE PATH "Path to MG-Tetra")
if(MG_TETRA)
    set(USE_MG_TETRA ON)
    message(STATUS "Using MG-Tetra directory =  ${MG_TETRA}")
endif(MG_TETRA)
 
# Optional custom configuration of data directories
if(RINGMESH_WITH_TESTS AND UNIX)
    list(APPEND CMAKE_CONFIGURATION_TYPES Coverage)
endif()

#------------------------------------------------------------------------------------------------
# Ensure consistent choice in configuration options
include(CMakeDependentOption)
# Option to build ringmesh-view. Disable when RINGMESH_WITH_GRAPHICS
# is OFF. The fist time RINGMESH_WITH_GRAPHICS is ON, this
# option is ON.
CMAKE_DEPENDENT_OPTION(BUILD_RINGMESH_VIEW "Build ringmesh-view" ON
    "RINGMESH_WITH_GRAPHICS" OFF)

# Check that RINGMesh has been cloned and not downloaded
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/.git)
    message(FATAL_ERROR 
        "This directory is not a complete git project." 
        "To configure RINGMesh, first clone it from its repository."
    )
endif()


# CMake is called, a main project (super build project) is generated first.
# This main project is responsible to create third party cmake projects including
# the one of RINGMesh. This is done during compilation, i.e., RINGMesh.cmake is
# called (USE_SUPERBUILD = OFF) when the super build project is compiled.
if(USE_SUPERBUILD)
    include(cmake/SuperBuild.cmake)
else(USE_SUPERBUILD)
    include(cmake/RINGMesh.cmake)
endif(USE_SUPERBUILD)

include(tools/external_tools.cmake)
