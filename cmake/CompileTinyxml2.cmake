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
# tinyxml2
# Set the path to tinyxml2 code
set(TINYXML2_PATH ${PROJECT_SOURCE_DIR}/third_party/tinyxml2)

# tinyxml2 platform dependent settings
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(TINYXML2_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/tinyxml2/${CMAKE_BUILD_TYPE})
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(TINYXML2_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/tinyxml2)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

# Define tinyxml2 as an external project that we know how to
# configure and compile
# Define Tinyxml2 as an external project that we know how to
# configure and compile
ExternalProject_Add(tinyxml2_ext
  PREFIX ${TINYXML2_PATH_BIN}

  #--Download step--------------
  DOWNLOAD_COMMAND ""

  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${TINYXML2_PATH}
      CONFIGURE_COMMAND ${CMAKE_COMMAND} ${TINYXML2_PATH}
          -G ${CMAKE_GENERATOR} 
          -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
          -DBUILD_TESTING:BOOL=OFF
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}

  #--Build step-----------------
  BINARY_DIR ${TINYXML2_PATH_BIN}
  #-- Command to build tinyxml2
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${TINYXML2_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_COMMAND ""
)

ExternalProject_Add_Step(tinyxml2_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add tinyxml2 include directories to the current ones
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party)

# Add tinyxml2 project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} debug tinyxml2d optimized tinyxml2)

# Add tinyxml2 bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${TINYXML2_PATH_BIN})
