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
# minizip
# Set the path to minizip code
set(MINIZIP_PATH ${PROJECT_SOURCE_DIR}/third_party/minizip)
set(MINIZIP_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/minizip/${CMAKE_BUILD_TYPE})
set(MINIZIP_INSTALL_PREFIX ${MINIZIP_PATH_BIN}/install/${CMAKE_CFG_INTDIR} CACHE INTERNAL "Minizip install directory")

# Define minizip as an external project that we know how to
# configure and compile
ExternalProject_Add(minizip_ext
  PREFIX ${MINIZIP_PATH_BIN}

  #--Download step--------------
  DOWNLOAD_COMMAND ""

  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${MINIZIP_PATH}
      CONFIGURE_COMMAND ${CMAKE_COMMAND} ${MINIZIP_PATH}
          -G ${CMAKE_GENERATOR}
          -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
          -DUSE_AES:BOOL=OFF
          -DZLIB_ROOT:PATH=${ZLIB_ROOT}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
          -DCMAKE_INSTALL_PREFIX:STRING=${MINIZIP_INSTALL_PREFIX}
          

  #--Build step-----------------
  BINARY_DIR ${MINIZIP_PATH_BIN}
  #-- Command to build minizip
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${MINIZIP_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_DIR ${MINIZIP_INSTALL_PREFIX}
  
  DEPENDS zlib_ext
)

ExternalProject_Add_Step(minizip_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add minizip include directories to the current ones
# same as minizip