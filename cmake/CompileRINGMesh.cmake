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
# RINGMesh
# Set the path to ringmesh code
set(RINGMesh_PATH ${PROJECT_SOURCE_DIR}/third_party/ringmesh)

# Define ringmesh as an external project that we know how to
# configure and compile
ExternalProject_Add(ringmesh_ext
  PREFIX ${PROJECT_BINARY_DIR}

  #--Download step--------------
  DOWNLOAD_COMMAND ""

  #--Update/Patch step----------
  UPDATE_COMMAND ""
  CMAKE_ARGS "-DZLIB_PATH_BIN=${ZLIB_PATH_BIN}"

  #--Configure step-------------
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
      CONFIGURE_COMMAND ${CMAKE_COMMAND} ${PROJECT_SOURCE_DIR}
          ${RINGMESH_EXTRA_ARGS} 
          -G ${CMAKE_GENERATOR}
          -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
          -DGLOBAL_BINARY_DIR=${PROJECT_BINARY_DIR}/..
          -DCMAKE_PROJECT_RINGMesh_INCLUDE=${PROJECT_SOURCE_DIR}/cmake/temp_fix.cmake
          -DUSE_SUPERBUILD=OFF

  #--Build step-----------------
  BINARY_DIR ${PROJECT_BINARY_DIR}
  #-- Command to build ringmesh
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_COMMAND ""
  
  DEPENDS geogram_ext tinyxml2_ext zlib_ext minizip_ext
)

ExternalProject_Add_Step(ringmesh_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

