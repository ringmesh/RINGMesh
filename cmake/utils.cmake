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

function(copy_for_windows target)

    # On windows, without proper installation steps, we need to
    # copy of Geogram dll and pdb information to RINGMesh
    # to be able to launch RINGMesh utilities and tests from the debugger

#    include(BundleUtilities)
    # The dll and debug info of RINGMesh are in
    # build/ringmesh/Debug or build/ringmesh/Release.
#    add_custom_command(TARGET ${target} POST_BUILD
#        COMMAND  fixup_bundle(\"$<TARGET_FILE_DIR:${target}>\" \"\" \"${GEOGRAM_PATH_BIN}/bin/$<CONFIGURATION>\")
#            COMMENT "Copy geogram binaries")
    #fixup_bundle(${target} "geogram" "${GEOGRAM_PATH_BIN}/bin/$<CONFIGURATION>")
#    get_target_property(directory ${target} RUNTIME_OUTPUT_DIRECTORY)
    #get_target_property(name ${target} LOCATION)
#    message(STATUS ${target})
#    message(STATUS "$<TARGET_FILE_NAME:${target}>")
if(WIN32)    
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${PROJECT_BINARY_DIR}/$<CONFIGURATION>"
            "$<TARGET_FILE_DIR:${target}>/$<CONFIGURATION>"
            COMMENT "Copy RINGMesh dll")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${GEOGRAM_PATH_BIN}/bin/$<CONFIGURATION>"
            "$<TARGET_FILE_DIR:${target}>/$<CONFIGURATION>"
            COMMENT "Copy geogram binaries")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${ZLIB_PATH_BIN}/$<CONFIGURATION>"
            "$<TARGET_FILE_DIR:${target}>/$<CONFIGURATION>"
            COMMENT "Copy zlib binaries")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${TINYXML2_PATH_BIN}/$<CONFIGURATION>"
            "$<TARGET_FILE_DIR:${target}>/$<CONFIGURATION>"
            COMMENT "Copy tinyxml2 binaries")
    add_custom_command(TARGET RINGMesh POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${MINIZIP_PATH_BIN}/$<CONFIGURATION>"
            "$<TARGET_FILE_DIR:${target}>/$<CONFIGURATION>"
            COMMENT "Copy minizip binaries")
endif(WIN32)
endfunction()

macro(add_ringmesh_executable exe_path folder_name)
    get_filename_component(exe_name ${exe_path} NAME_WE)

    # Set the target as an executable
    add_executable(${exe_name} ${exe_path})
    target_link_libraries(${exe_name} PRIVATE RINGMesh geogram)
    add_dependencies(${exe_name} RINGMesh)

    # Add the project to a folder of projects for the tests
    set_target_properties(${exe_name} PROPERTIES FOLDER ${folder_name})
    
    # ringmesh_files is defined in the root RINGMesh CMakeLists.txt.
    # This line is for clang utilities.
    set(ringmesh_files ${ringmesh_files} ${bin_path} PARENT_SCOPE)
    copy_for_windows(${exe_name})
endmacro()

function(add_ringmesh_binary bin_path)
    add_ringmesh_executable(${bin_path} "Utilities")
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endfunction()

function(add_ringmesh_utility bin_path)
    add_ringmesh_executable(${bin_path} "Utilities")
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/utilities)
endfunction()

function(add_ringmesh_test cpp_file_path)
    add_ringmesh_executable(${cpp_file_path} "Tests")
    # Add the test to CTest
    add_test(NAME ${exe_name} COMMAND ${exe_name})
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/tests)
endfunction()


function(add_ringmesh_tutorial cpp_file_path)
    add_ringmesh_test(${cpp_file_path} "Tutorials")
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/tutorials)
endfunction()
