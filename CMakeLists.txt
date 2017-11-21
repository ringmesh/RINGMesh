# Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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
set(USER_CONFIG "${CMAKE_SOURCE_DIR}/cmake/UserConfig.cmake" CACHE PATH
  "Path to optional user configuration file.")

#------------------------------------------------------------------------------------------------
# Options of the RINGMesh project
# They can be set
# 1 - using CMake GUI,
# 2 - in command line: cmake -D<var>:<type>=<value> MyProjectFolder
# 3 - by creating a UserConfig.cmake file in the cmake directory (see DefaultConfig.cmake)

# Load user settings before the defaults - this way the defaults will not
# overwrite the user set options. If the user has not set all options, we still
# have the defaults.
message(STATUS "(optional) USER_CONFIG = ${USER_CONFIG}")
include("${USER_CONFIG}" OPTIONAL)
include("${CMAKE_SOURCE_DIR}/cmake/DefaultConfig.cmake")

#------------------------------------------------------------------------------------------------
#ensure consistent choice in configuration options
include(CMakeDependentOption)
# Option to build ringmesh-view. Disable when RINGMESH_WITH_GRAPHICS
# is OFF. The fist time RINGMESH_WITH_GRAPHICS is ON, this
# option is ON.
CMAKE_DEPENDENT_OPTION(BUILD_RINGMESH_VIEW "Build ringmesh-view" ON
    "RINGMESH_WITH_GRAPHICS" OFF)

set(CMAKE_INSTALL_MESSAGE LAZY)

# CMake is called, a main project (super build project) is generated first.
# This main project is responsible to create third party cmake projects including
# the one of RINGMesh. This is done during compilation, i.e., RINGMesh.cmake is
# called (USE_SUPERBUILD = OFF) when the super build project is compiled.
if(USE_SUPERBUILD)
    include("${CMAKE_SOURCE_DIR}/cmake/SuperBuild.cmake")
else(USE_SUPERBUILD)
    include("${CMAKE_SOURCE_DIR}/cmake/RINGMesh.cmake")
endif(USE_SUPERBUILD)
