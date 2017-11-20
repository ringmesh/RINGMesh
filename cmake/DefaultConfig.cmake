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
#------------------------------------------------------------------------
# Default options. To replace these, copy this file to UserConfig.cmake and modify it.
#------------------------------------------------------------------------

# Configuration modes - Debug and Release by default + RelWithDebInfo for Windows.
# RelWithDebInfo is mandatory to debug some plugins in Windows in debug mode (e.g., a Gocad plugin) to
# avoid issues in particular in the STL (e.g., bad transfer of std::string from Gocad to RINGMesh).
set(CMAKE_CONFIGURATION_TYPES Debug Release RelWithDebInfo CACHE TYPE INTERNAL FORCE)

#------------------------------------------------------------------------------------------------
# Project options
# Optional components of RINGMesh - creation of specific targets
option(RINGMESH_WITH_GRAPHICS "Compile viewer" OFF)
option(RINGMESH_WITH_UTILITIES "Compile utility executables" OFF)
option(RINGMESH_WITH_TUTORIALS "Compile API trainings and tuturials" OFF)
option(RINGMESH_WITH_TESTS "Compile test projects" OFF)
# Optional custom steps 
option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)")
option(BUILD_GEOGRAM_WITHOUT_EXE "Compile geogram without tests and executables" ON)
mark_as_advanced(BUILD_GEOGRAM_WITHOUT_EXE)
option(USE_SUPERBUILD "Whether or not a superbuild should be invoked" ON)
mark_as_advanced(USE_SUPERBUILD)

#------------------------------------------------------------------------------------------------
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

