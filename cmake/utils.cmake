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

function(source_file_directory var directory)
    file(GLOB sources "${PROJECT_SOURCE_DIR}/src/ringmesh/${directory}/*.cpp")
    source_group("Source Files\\${directory}" FILES ${sources})
    set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

function(include_file_directory var directory)
    file(GLOB sources "${PROJECT_SOURCE_DIR}/include/ringmesh/${directory}/*.h")
    source_group("Header Files\\${directory}" FILES ${sources})
    set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

function(copy_for_windows)
if(WIN32)
    # On windows, without proper installation steps, we need to
    # copy of Geogram dll and pdb information to RINGMesh
    # to be able to launch RINGMesh utilities and tests from the debugger

    # The dll and debug info of RINGMesh are in
    # build/ringmesh/Debug or build/ringmesh/Release.
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${PROJECT_BINARY_DIR}/$<CONFIGURATION>"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>"
            COMMENT "Copy RINGMesh dll")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${GEOGRAM_PATH_BIN}/bin/$<CONFIGURATION>"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>"
            COMMENT "Copy geogram binaries")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${ZLIB_PATH_BIN}/$<CONFIGURATION>"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>"
            COMMENT "Copy zlib binaries")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${TINYXML2_PATH_BIN}/$<CONFIGURATION>"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>"
            COMMENT "Copy tinyxml2 binaries")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${MINIZIP_PATH_BIN}/$<CONFIGURATION>"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>"
            COMMENT "Copy minizip binaries")
endif(WIN32)
endfunction()

function(add_ringmesh_binary bin_path)
    get_filename_component(bin_name ${bin_path} NAME_WE)
    # Set the target as an executable
    add_executable(${bin_name} ${bin_path})

    set(BINARY_DEPENDANCIES RINGMesh geogram)
    foreach(arg ${ARGN})
        set(BINARY_DEPENDANCIES ${BINARY_DEPENDANCIES} ${arg})
    endforeach()
    target_link_libraries(${bin_name} ${BINARY_DEPENDANCIES})

    # Add the project to a folder of projects for the tests
    set_property(TARGET ${bin_name} PROPERTY FOLDER "Utilities")

    # ringmesh_files is defined in the root RINGMesh CMakeLists.txt.
    # This line is for clang utilities.
    set(ringmesh_files ${ringmesh_files} ${bin_path} PARENT_SCOPE)
endfunction()

function(add_ringmesh_test cpp_file_path folder_name)
    get_filename_component(cpp_file_name ${cpp_file_path} NAME_WE)

    # Add the target for that file
    add_executable(${cpp_file_name} ${cpp_file_path})

    # Link RINGMesh Make sure YourLib is linked to each app
    target_link_libraries(${cpp_file_name} geogram RINGMesh)

    add_dependencies(${cpp_file_name} RINGMesh)

    # Add the project to a folder of projects for the tests
    set_property(TARGET ${cpp_file_name} PROPERTY FOLDER ${folder_name})

    # Add the test to CTest
    add_test(NAME ${cpp_file_name} COMMAND ${cpp_file_name})
endfunction()
