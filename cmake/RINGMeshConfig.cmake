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

get_filename_component(RINGMESH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(RINGMESH_ROOT_DIRECTORY "${RINGMESH_CMAKE_DIR}/..")

set(RINGMesh_INCLUDE_DIR ${RINGMESH_ROOT_DIRECTORY}/include)
set(GEOGRAM_INCLUDE_DIR ${RINGMESH_ROOT_DIRECTORY}/third_party/geogram/src/lib)
set(THIRD_PARTY_INCLUDE_DIR ${RINGMESH_ROOT_DIRECTORY}/third_party)
set(THIRD_PARTY_INCLUDE_DIR 
	${THIRD_PARTY_INCLUDE_DIR} 
	${RINGMESH_ROOT_DIRECTORY}/third_party/zlib)
# On Visual Studio and Xcode there is no CMAKE_BUILD_TYPE variable since everything
# is in the same project. It is necessary to differenciate Release and Debug.
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # RINGMESH
    set(RINGMesh_CONFIG_INCLUDE_DIR ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/${CMAKE_BUILD_TYPE})
    find_library(RINGMesh_LIBRARY NAMES RINGMesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/${CMAKE_BUILD_TYPE}/lib)

    # GEOGRAM
    find_library(GEOGRAM_LIBRARY NAMES geogram
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/geogram/${CMAKE_BUILD_TYPE}/lib)
    find_library(GEOGRAM_GFX_LIBRARY NAMES geogram_gfx
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/geogram/${CMAKE_BUILD_TYPE}/lib)

    # ZLIB
    unset(ZLIB_LIBRARY CACHE)
    find_library(ZLIB_LIBRARY NAMES z
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/zlib/${CMAKE_BUILD_TYPE}
        NO_DEFAULT_PATH)

    # MINIZIP
    unset(MINIZIP_LIBRARY CACHE)
    find_library(MINIZIP_LIBRARY NAMES minizip minizipd
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/minizip/${CMAKE_BUILD_TYPE}
        NO_DEFAULT_PATH)

    # TINYXML2
    unset(TINYXML2_LIBRARY CACHE)
    find_library(TINYXML2_LIBRARY NAMES tinyxml2 tinyxml2d
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/tinyxml2/${CMAKE_BUILD_TYPE}
        NO_DEFAULT_PATH)
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # RINGMESH
    set(RINGMesh_CONFIG_INCLUDE_DIR ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh)
    find_library(RINGMesh_DEBUG_LIBRARY NAMES RINGMesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/lib/Debug)
    find_library(RINGMesh_RELEASE_LIBRARY NAMES RINGMesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/lib/Release)

    # GEOGRAM
    find_library(GEOGRAM_DEBUG_LIBRARY NAMES geogram
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/geogram/lib/Debug)
    find_library(GEOGRAM_GFX_DEBUG_LIBRARY NAMES geogram_gfx
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/geogram/lib/Debug)
    find_library(GEOGRAM_RELEASE_LIBRARY NAMES geogram
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/geogram/lib/Release)
    find_library(GEOGRAM_GFX_RELEASE_LIBRARY NAMES geogram_gfx
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/geogram/lib/Release)

    # ZLIB
    unset(ZLIB_DEBUG_LIBRARY CACHE)
    find_library(ZLIB_DEBUG_LIBRARY NAMES zlibd
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/zlib/Debug NO_DEFAULT_PATH)
    unset(ZLIB_RELEASE_LIBRARY CACHE)
    find_library(ZLIB_RELEASE_LIBRARY NAMES zlib
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/zlib/Release NO_DEFAULT_PATH)
    unset(MINIZIP_DEBUG_LIBRARY CACHE)

    # MINIZIP
    find_library(MINIZIP_DEBUG_LIBRARY NAMES minizipd
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/minizip/Debug NO_DEFAULT_PATH)
    unset(MINIZIP_RELEASE_LIBRARY CACHE)
    find_library(MINIZIP_RELEASE_LIBRARY NAMES minizip
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/minizip/Release NO_DEFAULT_PATH)
    unset(TINYXML2_DEBUG_LIBRARY CACHE)

    # TINYXML2
    find_library(TINYXML2_DEBUG_LIBRARY NAMES tinyxml2d
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/tinyxml2/Debug NO_DEFAULT_PATH)
    unset(TINYXML2_RELEASE_LIBRARY CACHE)
    find_library(TINYXML2_RELEASE_LIBRARY NAMES tinyxml2
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/third_party/tinyxml2/Release NO_DEFAULT_PATH)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

include(FindPackageHandleStandardArgs)
set(RINGMesh_INCLUDE_DIRS
    ${RINGMesh_INCLUDE_DIR}
    ${RINGMesh_CONFIG_INCLUDE_DIR}
    ${GEOGRAM_INCLUDE_DIR}
    ${THIRD_PARTY_INCLUDE_DIR})
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # Check that all the paths and libraries have been found.
    # Else error message in during cmake configuration.
    find_package_handle_standard_args(RINGMesh DEFAULT_MSG RINGMesh_INCLUDE_DIR
        RINGMesh_CONFIG_INCLUDE_DIR GEOGRAM_INCLUDE_DIR RINGMesh_LIBRARY
        GEOGRAM_LIBRARY ZLIB_LIBRARY TINYXML2_LIBRARY 
        MINIZIP_LIBRARY THIRD_PARTY_INCLUDE_DIR)
    set(RINGMesh_LIBRARIES ${RINGMesh_LIBRARY} ${GEOGRAM_LIBRARY} ${ZLIB_LIBRARY} )
    # Local variables are not displayed in the cmake-gui.
    # RINGMesh_DIR is defined by cmake implicitly.
    mark_as_advanced(RINGMesh_DIR RINGMesh_INCLUDE_DIR RINGMesh_CONFIG_INCLUDE_DIR
        GEOGRAM_INCLUDE_DIR THIRD_PARTY_INCLUDE_DIR RINGMesh_LIBRARY GEOGRAM_LIBRARY
         ZLIB_LIBRARY TINYXML2_LIBRARY MINIZIP_LIBRARY)
   if(NOT ${GEOGRAM_GFX_LIBRARY} STREQUAL "GEOGRAM_GFX_LIBRARY-NOTFOUND")
       set(RINGMesh_LIBRARIES ${RINGMesh_LIBRARIES} ${GEOGRAM_GFX_LIBRARY})
       mark_as_advanced(GEOGRAM_GFX_LIBRARY)
   endif(NOT ${GEOGRAM_GFX_LIBRARY} STREQUAL "GEOGRAM_GFX_LIBRARY-NOTFOUND")
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # Check that all the paths and libraries have been found.
    # Else error message in during cmake configuration.
    find_package_handle_standard_args(RINGMesh DEFAULT_MSG RINGMesh_INCLUDE_DIR
        RINGMesh_CONFIG_INCLUDE_DIR GEOGRAM_INCLUDE_DIR RINGMesh_RELEASE_LIBRARY
        RINGMesh_DEBUG_LIBRARY GEOGRAM_RELEASE_LIBRARY GEOGRAM_DEBUG_LIBRARY
        THIRD_PARTY_INCLUDE_DIR ZLIB_RELEASE_LIBRARY ZLIB_DEBUG_LIBRARY
        MINIZIP_RELEASE_LIBRARY MINIZIP_DEBUG_LIBRARY TINYXML2_RELEASE_LIBRARY
        TINYXML2_DEBUG_LIBRARY)
    # When the linking is done in RINGMesh dependent code, debug will just select
    # ${RINGMesh_DEBUG_LIBRARY} in RINGMesh_LIBRARIES for the Debug mode,
    # else (optimized) ${RINGMesh_RELEASE_LIBRARY} will be selected
    # (even for RelWithDebInfo...). The same for geogram.
    set(RINGMesh_LIBRARIES
        debug ${RINGMesh_DEBUG_LIBRARY}
        optimized ${RINGMesh_RELEASE_LIBRARY}
        debug ${GEOGRAM_DEBUG_LIBRARY}
        optimized ${GEOGRAM_RELEASE_LIBRARY}
        debug ${GEOGRAM_GFX_DEBUG_LIBRARY}
        optimized ${GEOGRAM_GFX_RELEASE_LIBRARY}
        debug ${ZLIB_DEBUG_LIBRARY}
        optimized ${ZLIB_RELEASE_LIBRARY} 
        debug ${MINIZIP_DEBUG_LIBRARY}
        optimized ${MINIZIP_RELEASE_LIBRARY} 
        debug ${TINYXML2_DEBUG_LIBRARY}
        optimized ${TINYXML2_RELEASE_LIBRARY} )
    # Local variables are not displayed in the cmake-gui.
    # RINGMesh_DIR is defined by cmake implicitly.
    mark_as_advanced(RINGMesh_DIR RINGMesh_INCLUDE_DIR RINGMesh_CONFIG_INCLUDE_DIR
        GEOGRAM_INCLUDE_DIR RINGMesh_RELEASE_LIBRARY RINGMesh_DEBUG_LIBRARY
        GEOGRAM_RELEASE_LIBRARY GEOGRAM_DEBUG_LIBRARY THIRD_PARTY_INCLUDE_DIR
        ZLIB_RELEASE_LIBRARY ZLIB_DEBUG_LIBRARY)
   if(NOT ${GEOGRAM_GFX_LIBRARY} STREQUAL "GEOGRAM_GFX_LIBRARY-NOTFOUND")
    set(RINGMesh_LIBRARIES
        ${RINGMESH_LIBRARIES}
        debug ${GEOGRAM_GFX_DEBUG_LIBRARY}
        optimized ${GEOGRAM_GFX_RELEASE_LIBRARY})
    mark_as_advanced(GEOGRAM_GFX_DEBUG_LIBRARY GEOGRAM_GFX_RELEASE_LIBRARY)
   endif(NOT ${GEOGRAM_GFX_LIBRARY} STREQUAL "GEOGRAM_GFX_LIBRARY-NOTFOUND")
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

function(set_ringmesh_includes_and_libs extra_libs)
if(RINGMesh_FOUND)
    message(STATUS "RINGMesh is found!")
    message(STATUS "RINGMesh include directory ${RINGMesh_INCLUDE_DIR}")
    message(STATUS "geogram include directory ${GEOGRAM_INCLUDE_DIR}")
    message(STATUS "third party include directory ${THIRD_PARTY_INCLUDE_DIR}")
    message(STATUS "RINGMesh include directories ${RINGMesh_INCLUDE_DIRS}")
    if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
        message(STATUS "RINGMesh library ${RINGMesh_LIBRARY}")
        message(STATUS "geogram library ${GEOGRAM_LIBRARY}")
        message(STATUS "zlib library ${ZLIB_LIBRARY}")
    else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
        message(STATUS "RINGMesh release library ${RINGMesh_RELEASE_LIBRARY}")
        message(STATUS "geogram release library ${GEOGRAM_RELEASE_LIBRARY}")
        message(STATUS "zlib release library ${ZLIB_RELEASE_LIBRARY}")
        message(STATUS "RINGMesh debug library ${RINGMesh_DEBUG_LIBRARY}")
        message(STATUS "geogram debug library ${GEOGRAM_DEBUG_LIBRARY}")
        message(STATUS "zlib debug library ${ZLIB_DEBUG_LIBRARY}")
    endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    message(STATUS "RINGMesh libraries ${RINGMesh_LIBRARIES}")
	
    # Add RINGMesh include directories
    include_directories(${RINGMesh_INCLUDE_DIRS})
    # Add RINGMesh project libs to the libs with which the client code will link
    set(${extra_libs} ${${extra_libs}} ${RINGMesh_LIBRARIES} PARENT_SCOPE)
else(RINGMesh_FOUND)
    # In theory find_package with REQUIRED stops the cmake generation if the package is not found, so this else should never happen...
    message(FATAL_ERROR "RINGMesh not found")
endif(RINGMesh_FOUND)
endfunction()
