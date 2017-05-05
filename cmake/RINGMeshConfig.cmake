# Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
# All rights reserved.
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
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#     http://www.ring-team.org
# 
#     RING Project
#     Ecole Nationale Superieure de Geologie - GeoRessources
#     2 Rue du Doyen Marcel Roubault - TSA 70605
#     54518 VANDOEUVRE-LES-NANCY
#     FRANCE

get_filename_component(RINGMESH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
SET(RINGMESH_ROOT_DIRECTORY "${RINGMESH_CMAKE_DIR}/..")

find_path(RINGMesh_INCLUDE_DIR NAMES ringmesh
    PATHS ${RINGMESH_ROOT_DIRECTORY}/include)
# On Windows there is no CMAKE_BUILD_TYPE variable since everything
# is in the same project. It is necessary to differenciate Release and Debug.
if(WIN32)
    find_path(RINGMesh_CONFIG_INCLUDE_DIR NAMES ringmesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh)
    find_library(RINGMesh_DEBUG_LIBRARY NAMES RINGMesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/lib/Debug)
    find_library(RINGMesh_RELEASE_LIBRARY NAMES RINGMesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/lib/Release)
    find_library(GEOGRAM_DEBUG_LIBRARY NAMES geogram
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/geogram/lib/Debug)
    find_library(GEOGRAM_RELEASE_LIBRARY NAMES geogram
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/geogram/lib/Release)
else(WIN32)
    find_path(RINGMesh_CONFIG_INCLUDE_DIR NAMES ringmesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/${CMAKE_BUILD_TYPE})
    find_library(RINGMesh_LIBRARY NAMES RINGMesh
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/ringmesh/${CMAKE_BUILD_TYPE}/lib)
    find_library(GEOGRAM_LIBRARY NAMES geogram
        PATHS ${RINGMESH_ROOT_DIRECTORY}/build/geogram/${CMAKE_BUILD_TYPE}/lib)
endif(WIN32)
find_path(GEOGRAM_INCLUDE_DIR NAMES geogram
    PATHS ${RINGMESH_ROOT_DIRECTORY}/third_party/geogram/src/lib)
find_path(THIRD_PARTY_INCLUDE_DIR NAMES zlib
    PATHS ${RINGMESH_ROOT_DIRECTORY}/third_party)


include(FindPackageHandleStandardArgs)
set(RINGMesh_INCLUDE_DIRS
    ${RINGMesh_INCLUDE_DIR}
    ${RINGMesh_CONFIG_INCLUDE_DIR}
    ${GEOGRAM_INCLUDE_DIR}
    ${THIRD_PARTY_INCLUDE_DIR})
if(WIN32)
    # Check that all the paths and libraries have been found.
    # Else error message in during cmake configuration.
    find_package_handle_standard_args(RINGMesh DEFAULT_MSG RINGMesh_INCLUDE_DIR
        RINGMesh_CONFIG_INCLUDE_DIR GEOGRAM_INCLUDE_DIR RINGMesh_RELEASE_LIBRARY
        RINGMesh_DEBUG_LIBRARY GEOGRAM_RELEASE_LIBRARY GEOGRAM_DEBUG_LIBRARY
        THIRD_PARTY_INCLUDE_DIR)
    # When the linking is done in RINGMesh dependent code, debug will just select
    # ${RINGMesh_DEBUG_LIBRARY} in RINGMesh_LIBRARIES for the Debug mode,
    # else (optimized) ${RINGMesh_RELEASE_LIBRARY} will be selected
    # (even for RelWithDebInfo...). The same for geogram.
    set(RINGMesh_LIBRARIES
        debug ${RINGMesh_DEBUG_LIBRARY}
        optimized ${RINGMesh_RELEASE_LIBRARY}
        debug ${GEOGRAM_DEBUG_LIBRARY}
        optimized ${GEOGRAM_RELEASE_LIBRARY})
    # Local variables are not displayed in the cmake-gui.
    # RINGMesh_DIR is defined by cmake implicitly.
    mark_as_advanced(RINGMesh_DIR RINGMesh_INCLUDE_DIR RINGMesh_CONFIG_INCLUDE_DIR
        GEOGRAM_INCLUDE_DIR RINGMesh_RELEASE_LIBRARY RINGMesh_DEBUG_LIBRARY
        GEOGRAM_RELEASE_LIBRARY GEOGRAM_DEBUG_LIBRARY THIRD_PARTY_INCLUDE_DIR)
else(WIN32)
    # Check that all the paths and libraries have been found.
    # Else error message in during cmake configuration.
    find_package_handle_standard_args(RINGMesh DEFAULT_MSG RINGMesh_INCLUDE_DIR
        RINGMesh_CONFIG_INCLUDE_DIR GEOGRAM_INCLUDE_DIR RINGMesh_LIBRARY
        GEOGRAM_LIBRARY THIRD_PARTY_INCLUDE_DIR)
    set(RINGMesh_LIBRARIES ${RINGMesh_LIBRARY} ${GEOGRAM_LIBRARY} )
    # Local variables are not displayed in the cmake-gui.
    # RINGMesh_DIR is defined by cmake implicitly.
    mark_as_advanced(RINGMesh_DIR RINGMesh_INCLUDE_DIR RINGMesh_CONFIG_INCLUDE_DIR
        GEOGRAM_INCLUDE_DIR THIRD_PARTY_INCLUDE_DIR RINGMesh_LIBRARY GEOGRAM_LIBRARY)
endif(WIN32)

