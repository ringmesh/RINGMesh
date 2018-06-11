# Copyright (c) 2018, Association Scientifique pour la Geologie et ses
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

SET(FESAPI_PATH ${PROJECT_SOURCE_DIR}/third_party/fesapi)
SET(FESAPI_PATH_BIN ${PROJECT_BINARY_DIR}/third_party/fesapi)
SET(FESAPI_INSTALL_PREFIX ${FESAPI_PATH_BIN}/install
    CACHE INTERNAL "Fesapi install directory")

SET (HDF5_ZIP_DESTIN_DIR "${PROJECT_BINARY_DIR}/third_party/hdf5/")
SET (HDF5_ZIP_DESTIN "${HDF5_ZIP_DESTIN_DIR}/hdf5.tar.gz")
SET (HDF5_ZIP_ROOT "hdf5-1.8.20-linux-centos7-x86_64-gcc485-shared")

IF(NOT EXISTS ${HDF5_ZIP_DESTIN})
    SET (HDF5_ZIP "${HDF5_ZIP_ROOT}.tar.gz")
    SET (HDF5_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/bin/${HDF5_ZIP}")

    MESSAGE(STATUS "Downloading: ${HDF5_URL}")
    FILE(MAKE_DIRECTORY ${HDF5_ZIP_DESTIN_DIR})
    FILE(DOWNLOAD ${HDF5_URL} ${HDF5_ZIP_DESTIN})
    MESSAGE(STATUS "Downloaded: ${HDF5_URL}")

    EXECUTE_PROCESS(
        COMMAND tar -C ${HDF5_ZIP_DESTIN_DIR} -xzf ${HDF5_ZIP_DESTIN}
    )
ENDIF()

SET (HDF5_INSTALL_PREFIX "${HDF5_ZIP_DESTIN_DIR}/${HDF5_ZIP_ROOT}/"
    CACHE PATH "Path to the directory which contains HDF5")
MARK_AS_ADVANCED(HDF5_INSTALL_PREFIX)

#MINIZIP_1_1
SET (MINIZIP_1_1_DESTIN_DIR "${PROJECT_BINARY_DIR}/third_party/minizip_1_1/")
SET (MINIZIP_1_1_DESTIN "${MINIZIP_1_1_DESTIN_DIR}/minizip_1_1.tar.gz")
SET (MINIZIP_1_1_ROOT "minizip-1.1-linux-ubuntu1604-x86_64-gcc540")

IF(NOT EXISTS ${MINIZIP_1_1_DESTIN})
    SET (MINIZIP_1_1 "${MINIZIP_1_1_ROOT}.tar.gz")
    SET (MINIZIP_1_1_URL "https://github.com/F2I-Consulting/Minizip/releases/download/1.1/${MINIZIP_1_1}")

    MESSAGE(STATUS "Downloading: ${MINIZIP_1_1_URL}")
    FILE(MAKE_DIRECTORY ${MINIZIP_1_1_DESTIN_DIR})
    FILE(DOWNLOAD ${MINIZIP_1_1_URL} ${MINIZIP_1_1_DESTIN})
    MESSAGE(STATUS "Downloaded: ${MINIZIP_1_1_URL}")

    EXECUTE_PROCESS(
        COMMAND tar -C ${MINIZIP_1_1_DESTIN_DIR} -xzf ${MINIZIP_1_1_DESTIN}
    )
ENDIF()

SET (MINIZIP_1_1_INSTALL_PREFIX "${MINIZIP_1_1_DESTIN_DIR}/${MINIZIP_1_1_ROOT}/"
    CACHE PATH "Path to the directory which contains MINIZIP_1_1")
MARK_AS_ADVANCED(MINIZIP_1_1_INSTALL_PREFIX)

ExternalProject_Add(fesapi_ext
    PREFIX ${FESAPI_PATH_BIN}
    SOURCE_DIR ${FESAPI_PATH}
    CMAKE_GENERATOR ${CMAKE_GENERATOR}
    CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_MESSAGE=LAZY
    CMAKE_CACHE_ARGS
        -DBUILD_TESTING:BOOL=OFF
        -DCMAKE_INSTALL_PREFIX:PATH=${FESAPI_INSTALL_PREFIX}
        -DCMAKE_INSTALL_LIBDIR:PATH=lib
        -DCMAKE_MACOSX_RPATH:BOOL=ON
        -DZLIB_ROOT:PATH=${ZLIB_ROOT}
        -DHDF5_INSTALL_PREFIX:PATH=${HDF5_INSTALL_PREFIX}
        -DMINIZIP_1_1_INSTALL_PREFIX:PATH=${MINIZIP_1_1_INSTALL_PREFIX}
    BINARY_DIR ${FESAPI_PATH_BIN}
    BUILD_ALWAYS 1
    INSTALL_DIR ${FESAPI_INSTALL_PREFIX}
    DEPENDS zlib_ext
)
