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

set(GEOGRAM_PATH ${PROJECT_SOURCE_DIR}/third_party/geogram)
set(GEOGRAM_PATH_BIN ${PROJECT_BINARY_DIR}/third_party/geogram)
set(GEOGRAM_INSTALL_PREFIX ${GEOGRAM_PATH_BIN}/install
    CACHE INTERNAL "Geogram install directory")

if(WIN32)
    set(geoplatform Win-vs-dynamic-generic)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        if(APPLE)
            set(geoplatform Darwin-clang-dynamic)
        else()
            set(geoplatform Linux64-clang-dynamic)
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(geoplatform Linux64-gcc-dynamic)
    endif()
endif()

ExternalProject_Add(geogram_ext
    PREFIX ${GEOGRAM_PATH_BIN}
    SOURCE_DIR ${GEOGRAM_PATH}
    CMAKE_GENERATOR ${CMAKE_GENERATOR}
    CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_MESSAGE=LAZY
    CMAKE_CACHE_ARGS
        -DVORPALINE_PLATFORM:STRING=${geoplatform}
        -DGEOGRAM_WITH_LUA:BOOL=OFF
        -DGEOGRAM_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN} 
        -DGEOGRAM_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
        -DGEOGRAM_WITH_EXPLORAGRAM:BOOL=OFF
        -DGEOGRAM_LIB_ONLY:BOOL=${BUILD_GEOGRAM_WITHOUT_EXE}
        -DCMAKE_INSTALL_PREFIX:PATH=${GEOGRAM_INSTALL_PREFIX}
    BINARY_DIR ${GEOGRAM_PATH_BIN}
    BUILD_ALWAYS 1
    INSTALL_DIR ${GEOGRAM_INSTALL_PREFIX}
)
