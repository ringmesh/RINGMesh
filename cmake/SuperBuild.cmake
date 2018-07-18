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
# Get all the submodules
set(submodules data third_party/zlib third_party/tinyxml2)
if(RINGMESH_WITH_GUI)
    set(submodules ${submodules} third_party/nbind)
endif()
execute_process(
   COMMAND git submodule update --init ${submodules}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set(RINGMESH_EXTRA_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    -DMG_TETRA:STRING=${MG_TETRA}
    -DRINGMESH_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN}
    -DRINGMESH_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
    -DBUILD_RINGMESH_VIEW:BOOL=${BUILD_RINGMESH_VIEW}
    -DRINGMESH_WITH_UTILITIES:BOOL=${RINGMESH_WITH_UTILITIES}
    -DRINGMESH_WITH_TESTS:BOOL=${RINGMESH_WITH_TESTS}
    -DRINGMESH_WITH_TUTORIALS:BOOL=${RINGMESH_WITH_TUTORIALS}
    -DBUILD_GEOGRAM_WITHOUT_EXE:BOOL=${BUILD_GEOGRAM_WITHOUT_EXE}
    -DRINGMESH_WITH_GUI:BOOL=${RINGMESH_WITH_GUI}
    -DRINGMESH_WITH_METIS:BOOL=${RINGMESH_WITH_METIS}
)

if(CPACK_PACKAGE_FILE_NAME)
    set(RINGMESH_EXTRA_ARGS 
        ${RINGMESH_EXTRA_ARGS}
        -DCPACK_PACKAGE_FILE_NAME:STRING=${CPACK_PACKAGE_FILE_NAME}
    )
endif()

#------------------------------------------------------------------------------------------------
# Generate configuration directories for single-configuration generators (Make)
# and run cmake configuration command in each one of them
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # If the CMAKE_BUILD_TYPE is no defined
    # i.e. at first run of that file
    if(NOT CMAKE_BUILD_TYPE)
      # For all configuration types
      foreach(config IN LISTS CMAKE_CONFIGURATION_TYPES)
          project(MULTI_CONFIG NONE)
          # The binary directory for this configuration
          set(project_binary_dir_config  ${PROJECT_BINARY_DIR}/${config})

          # Create the directory
          file(MAKE_DIRECTORY ${project_binary_dir_config})

          # Launch cmake for this configuration
          # by specifying the CMAKE_BUILD_TYPE and the PROJECT _BINARY_DIR,
          # and pass on all the options
          execute_process(
             COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}
                -G ${CMAKE_GENERATOR}
                ${RINGMESH_EXTRA_ARGS} 
                -DCMAKE_BUILD_TYPE=${config}
                WORKING_DIRECTORY ${project_binary_dir_config})
      endforeach()

      # Get out when all configurations have been run
      # We do not want to execute the rest of the file in that case
      return()
    endif()
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

# Execute the superbuild
project(SUPERBUILD)

# Additional cmake modules
include(ExternalProject)

include(${PROJECT_SOURCE_DIR}/cmake/CompileGeogram.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileTinyxml2.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileZlib.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileMinizip.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileRINGMesh.cmake)
